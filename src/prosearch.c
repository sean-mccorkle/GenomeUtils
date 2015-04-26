/* Program:     prosearch.c                                                  */
/* Programmer:  Sean R. McCorkle                                             */
/*              Biology Dept., Brookhaven National Laboratory                */
/* Language:    C                                                            */
/*                                                                           */
/* Description: Search DNA sequence files (FASTA format) for candidate       */
/*              promoter sequences, with adjustable mismatch threshold       */
/*                                                                           */
/* Usage:       prosearch [-aBmNSFvhV] <pat> [<file> [...]]                  */
/*                                                                           */
/*              where <file>s are DNA sequence files (FASTA format).  If no  */
/*              files are given, stdin is scanned.  "-" may also be used as  */
/*              a file name to indicate stdin.                               */
/*                                                                           */
/*              <pat>  is a short DNA sequence composed of any of            */
/*                                                                           */
/*                uppercase A,C,G,T - which will always match the appropriate*/
/*                                    position - they will never mismatch    */
/*                lowercase a,c,g,t - which match the appropriate base, but  */
/*                                    which are allowed to mismatch provided */
/*                                    the total number of mismatches doesn't */
/*                                    exceed the value specifed by -m        */
/*                ., n, N (wildcard)- matches any base and never counts as a */
/*                                    mismatch.                              */
/*                M,R,W,S,Y,K,V,H,D,B- match appropriate bases ONLY and are  */
/*                                    not allowed to mismatch outside their  */
/*                                    definitions                            */
/*                m,r,w,s,y,k,v,h,d,b- match appropriate bases and ARE       */
/*                                    allowed to mismatch outside their      */
/*                                    definitions, provided that the total   */
/*                                    number of mismatches doesn't exceed the*/
/*                                    value specifed by -m                   */
/*                                                                           */
/* Options:     -a<n>  accept up to <n> ambiguity codes in the search        */
/*                     sequence for any match (either upper or lower case    */
/*                     in the pattern is matched)                            */
/*              -B<n>  Bisulfite modify the DNA before matching              */
/*                      -B1    change all C -> T except C's in CpGs          */
/*                      -B2    change all C -> T including C's in CpGs       */
/*              -m<n>  accept up to <n> mismatches in LOWER case bases in    */
/*                     <pat>.  Default is 0 (exact matches).                 */
/*              -N<n>   print neighborhood <n> on each side                  */
/*              -S     print sequence names                                  */
/*              -F     print filenames                                       */
/*              -v     verbose output                                        */
/*              -h     print help, then exit                                 */
/*              -V     print version, then exit                              */
/*                                                                           */
/* Example:       prosearch -m1  GCAcct.ac                                   */
/*                                                                           */
/*                   will match GCACCTAAC                                    */
/*                              GCACCTCAC                                    */
/*                              GCACCTGAC                                    */
/*                              GCACCTTAC                                    */
/*                   and also will match any of these with one additional    */
/*                   mismatch that doesn't occur in the first three positions*/
/*                                                                           */
/* Caveats:                                                                  */
/*             After v1.3, ambiguity codes in the input sequence are         */
/*             accepted, up to the limit set by -a.  However, only the       */
/*             FIRST possible match is reported.  This may impact any        */
/*             statistical analysis performed on the output of this program  */
/*                                                                           */
/*                                                                           */
/*****************************************************************************/

#include <ctype.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>


#define MAX_HDR_LEN      256
#define MAX_PAT_LEN      256
#define MAX_STR_LEN      256
#define MAX_NEIGHBOR_LEN 50
#define MAX_BUFF_LEN     (MAX_NEIGHBOR_LEN  + MAX_PAT_LEN + MAX_NEIGHBOR_LEN)

#ifndef NULL
#define NULL        ((void *) 0)
#endif

char  bases[] = "ACGT";
/* char  templ[] = "ttcaGCacc.cGGacagc.cc"; */
char  templ[MAX_PAT_LEN+1];
int   templ_len;

int   buff_len;

int   amb_thresh      = 0;   /* set by -a */
int   bisulfite_level = 0;   /* set by -B */
int   mis_thresh      = 0;   /* set by -m */
int   neighbor_len    = 0;   /* set by -N */
int   print_headers   = 0;   /* set by -S */
int   print_filenames = 0;   /* set by -F */
int   verbose         = 0;   /* set by -v */

#ifdef SHOW_PERM_STATS
size_t       num_perms = 0;
#endif

void  version( void )                      /* print version number */
   {
    char *v;
    v = strdup("$Revision: 1.9 $"+11 );
    v[strlen(v)-2] = '\0';
    printf( "prosearch v%s\n", v );
   }


void  help( void )
   {
    version();
    printf( " \n\
             Search DNA sequence files (FASTA format) for candidate         \n\
             promoter sequences, with adjustable mismatch threshold.        \n\
                                                                            \n\
Usage:       prosearch [-aBmNSFvhV]  <pat> [<file> [...]]                   \n\
                                                                            \n\
             where <file>s are DNA sequence files (FASTA format).  If no    \n\
             files are given, stdin is scanned.  \"-\" may also be used as  \n\
             a file name to indicate stdin.                                 \n\
                                                                            \n\
             <pat>  is a short DNA sequence composed of any of              \n\
                                                                            \n\
               uppercase A,C,G,T   - which will always match the appropriate\n\
                                     position - they will never mismatch    \n\
               lowercase a,c,g,t   - which match the appropriate base, but  \n\
                                     which are allowed to mismatch provided \n\
                                     the total number of mismatches doesn't \n\
                                     exceed the value specifed by -m        \n\
               ., n, N (wildcard)  - matches any base and never counts as a \n\
                                     mismatch.                              \n\
               M,R,W,S,Y,K,V,H,D,B - match appropriate bases ONLY and are   \n\
                                     not allowed to mismatch outside their  \n\
                                     definitions                            \n\
               m,r,w,s,y,k,v,h,d,b - match appropriate bases and ARE        \n\
                                     allowed to mismatch outside their      \n\
                                     definitions, provided that the total   \n\
                                     number of mismatches doesn't exceed the\n\
                                     value specifed by -m                   \n\
                                                                            \n\
Options:     -a<n>  accept up to <n> ambiguity codes in the search          \n\
                    sequence for any match (either upper or lower case      \n\
                    in the pattern is matched)                              \n\
             -B<n>  Bisulfite modify the DNA before matching                \n\
                     -B1    change all C -> T except C's in CpGs            \n\
                     -B2    change all C -> T including C's in CpGs         \n\
             -m<n>  accept up to <n> mismatches in LOWER case bases in      \n\
                    <pat>.  Default is 0 (exact matches).                   \n\
             -N<n>  print neighboring sequences of length <n>               \n\
             -S     print sequence names                                    \n\
             -F     print filenames                                         \n\
             -v     verbose output                                          \n\
             -h     print help, then exit                                   \n\
             -V     print version, then exit                                \n\
                                                                            \n\
Example:       prosearch -m1  GCAcct.ac                                     \n\
                                                                            \n\
                  will match GCACCTAAC                                      \n\
                             GCACCTCAC                                      \n\
                             GCACCTGAC                                      \n\
                             GCACCTTAC                                      \n\
                  and also will match any of these with one additional      \n\
                  mismatch that doesn't occur in the first three positions  \n\
                                                                            \n\
Caveats:                                                                    \n\
              Ambiguity codes in the input sequence files are NOT handled   \n\
              properly.  They are ignored (for now).                        \n\
\n\
");
   } 

/* verfies that val is inside [a,b] inclusive.  Dies with an error  */
/* message if not                                                   */

void  check_int_range( int val, int a, int b, char *s )
   {
    if ( val < a || val > b )
       {
        fprintf( stderr, "%s must be between %d and %d (inclusive)\n",
                           s, a, b );
        fprintf( stderr, "run \"proserearch -h\" for more help\n" );
        exit( 1 );
       }
   }


void  parse_args( int argc, char **argv, int *nfiles, char ***filenames )
   {
    extern char *optarg;
    extern int   optind;
    int          c;
    static char *def_files[] = { "-", "" };

    while ( (c = getopt( argc, argv, "a:B:hm:FN:SvV" ) ) != -1 )
        switch( c )
           {
            case  'a':   amb_thresh = atoi( optarg );
                         check_int_range( amb_thresh, 1, MAX_PAT_LEN, 
                                          "-a value" );
                         break;
            case  'B':   bisulfite_level = atoi( optarg );
                         check_int_range( bisulfite_level, 0, 2, "-B value" );
                         break;
            case  'h':   help();
                         exit( 0 );
            case  'm':   mis_thresh = atoi( optarg );
                         check_int_range( mis_thresh, 1, MAX_PAT_LEN, 
                                          "-m value" );
                         break;
            case  'N':   neighbor_len = atoi( optarg );
                         check_int_range( neighbor_len, 1, MAX_NEIGHBOR_LEN, 
                                         "-N value" );
                         break; 
            case  'S':   print_headers = 1;
                         break;
            case  'F':   print_filenames = 1;
                         break;
            case  'v':   verbose = 1;
                         break;
            case  'V':   version();
                         exit( 0 );
            default:     help();
                         exit( 1 );
           }
    argc -= optind;
    argv += optind;
    if ( argc <= 0 )
       {
        fprintf( stderr, "no pattern specified.  (prosearch -h for help)\n" );
        exit( 1 );
       }
    /* Note to self:  needed here is a length check and a check on allowable */
    /* characters in pattern                                                 */
    strncpy( templ, *argv, MAX_PAT_LEN );
    argc--;
    argv++;
    if ( argc > 0 )
       {
        *nfiles = argc;
        *filenames = argv;
       }
    else
       {
        *nfiles = 1;
        *filenames = def_files;
       }
   }

/* comp_char[c] contains the complementary nucleotide character for c */

char comp_char[] = {
                  /*          0     1     2     3     4     5     6     7   */
                  /* 000 */  'X',  'X',  'X',  'X',  'X',  'X',  'X',  'X',
                  /* 010 */  'X',  'X',  'X',  'X',  'X',  'X',  'X',  'X',
                  /* 020 */  'X',  'X',  'X',  'X',  'X',  'X',  'X',  'X',
                  /* 030 */  'X',  'X',  'X',  'X',  'X',  'X',  'X',  'X',
                  /* 040 */  'X',  'X',  'X',  'X',  'X',  'X',  'X',  'X',
                  /* 050 */  'X',  'X',  'X',  'X',  'X',  'X',  'X',  'X',
                  /* 060 */  'X',  'X',  'X',  'X',  'X',  'X',  'X',  'X',
                  /* 070 */  'X',  'X',  'X',  'X',  'X',  'X',  'X',  'X',
                  /* 100 */  'X',  'T',  'V',  'G',  'H',  'X',  'X',  'C',
                  /* 110 */  'D',  'X',  'X',  'M',  'X',  'K',  'N',  'X',
                  /* 120 */  'X',  'X',  'Y',  'S',  'A',  'X',  'B',  'W',
                  /* 130 */  'X',  'R',  'X',  'X',  'X',  'X',  'X',  'X',
                  /* 140 */  'X',  't',  'v',  'g',  'h',  'X',  'X',  'c',
                  /* 150 */  'd',  'X',  'X',  'm',  'X',  'k',  'n',  'X',
                  /* 160 */  'X',  'X',  'y',  's',  'a',  'X',  'b',  'w',
                  /* 170 */  'X',  'r',  'X',  'X',  'X',  'X',  'X',  'X',
                 };

void  rc( char *r, char *s )
   {
    char *p;

    p = s + strlen( s );
    while ( p > s )
       *r++ = comp_char[*(--p)];
    *r = '\0';
   }


/* allowed_matches maps each nucleotide and ambiguity code to an array  */
/* of all the nucleotides (ACGT) that will match it.                    */

char *allowed_matches[256]; 

/* this initializes the allowed_matches table */

void  fill_allowed_matches( void )
   {
    int i;

    for ( i = 0; i < 256; i++ )      /* all characters map to NULL */
        allowed_matches[i] = NULL;   /* except for the DNA codes: */

    allowed_matches['a'] = allowed_matches['A'] = "A";
    allowed_matches['c'] = allowed_matches['C'] = "C";
    allowed_matches['g'] = allowed_matches['G'] = "G";
    allowed_matches['t'] = allowed_matches['T'] = "T";

    allowed_matches['m'] = allowed_matches['M'] = "AC";
    allowed_matches['r'] = allowed_matches['R'] = "AG";
    allowed_matches['w'] = allowed_matches['W'] = "AT";
    allowed_matches['s'] = allowed_matches['S'] = "CG";
    allowed_matches['y'] = allowed_matches['Y'] = "CT";
    allowed_matches['k'] = allowed_matches['K'] = "GT";

    allowed_matches['v'] = allowed_matches['V'] = "ACG";
    allowed_matches['h'] = allowed_matches['H'] = "ACT";
    allowed_matches['d'] = allowed_matches['D'] = "AGT";
    allowed_matches['b'] = allowed_matches['B'] = "GCT";

    allowed_matches['n'] = allowed_matches['N'] = "ACGT";
   }


/**************************************************************************/
/* match( a, b ) compares a, which must be a nucleotide A, C, G or T with */
/* b, which may be a nucleotide or an ambiguity code.  If a is in one     */
/* of the allowed matches of b as listed above, 1 is returned, otherwise  */
/* 0 is returned, indicating a mismatch                                   */
/**************************************************************************/

int  match( char a, char b )
   {
    char c;
    char *matches;
    int   k, n;

    c = toupper( a );
    matches = allowed_matches[b];
    n = strlen( matches );
    for ( k = 0; k < n; k++ )
        if ( c == matches[k] )
            return( 1 );
    return( 0 );
   }


                          /********************/
                          /* Permutation tree */
                          /********************/

int ind[128];  /* maps a,g,c,t chars to 0,1,2,3 indices, set in init_tree() */

char nt[] = { 'A', 'C', 'G', 'T' };  /* reverse map of ind[] */

typedef struct tnode {
                       struct tnode *ch[4];
                       char         *desc;
                     } TNODE;

TNODE *t_root;            /* tree root */

TNODE *new_tnode( void )
   {
    TNODE *t;

    if ( !( t = (TNODE *) malloc( sizeof(TNODE) ) ) )
       {
        perror( "can't allocate a tree node" );
        exit( errno );
       }
    t->ch[0] = t->ch[1] = t->ch[2] = t->ch[3] = NULL;
    t->desc = NULL;
    return( t );
   }

/*******************/
/* initialize tree */
/*******************/

void  init_tree( void )
   {
    int i;

    t_root = new_tnode();

    for ( i = 0; i < 128; i++ )
        ind[i] = -1;   
    ind['a'] = ind['A'] = 0;
    ind['c'] = ind['C'] = 1;
    ind['g'] = ind['G'] = 2;
    ind['t'] = ind['T'] = 3;
    if ( verbose )
       printf( "tree created\n" );
   }


/************************************************************************/
/* search tree for string.  If found, return value string for the match */
/* match, otherwise return NULL                                         */
/************************************************************************/

char *tree_lookup( char *s, int len )
   {
    TNODE *t;

    if ( verbose )
        printf( "tree_lookup [%s]\n", s );
    t = t_root;
    while ( len > 0 && ind[*s] >= 0 && t->ch[ind[*s]] !=NULL )
       {
        t = t->ch[ind[*s++]];
        len--;
       }

    if ( len == 0 )
        return( t->desc );
    else if ( ind[*s] < 0 )
        return( NULL );       /* 'X' or anything else in s fails */
    else
        return( NULL );
   }


char  *tree_rlookup( TNODE *t, char *s, int len, int allowed_ambs )
   {
    TNODE  *r;
    int     i, k, n;
    char   *matches;
    char   *res;
    
    if ( len <= 0 )                           /* end of string, so we're done*/
        return( t->desc );                    /* return the result           */
    else if ( (i = ind[*s]) >= 0 )            /* is this a nucl.? (a,c,g,t)? */
       {                                      /* if so,  is this nucleotide  */
        if ( (r = t->ch[i]) != NULL )         /* present here in the tree?   */
                                              /* then recurse...             */
            return( tree_rlookup( r, s+1, len-1, allowed_ambs ) );
        else
            return( NULL );                   /* otherwise stop: not found */
       }
    else if ( (matches = allowed_matches[*s]) ) /* is this an ambiguity code? */
       {
        if ( allowed_ambs > 0 )               /* can we match one more amb? */
           {
            n = strlen( matches );            /* then for each possibility, */
            for ( k = 0; k < n; k++ )         /* if its in the tree, recurse*/
                if ( ( r = t->ch[ind[matches[k]]] ) ) /* and return 1st hit */
                    if ( (res = tree_rlookup( r, s+1, len-1, allowed_ambs-1 ) ) )
                        return( res );
            return( NULL );                   /* not found if we finished loop*/
           }
        else                                  /* no more ambiguities allowed, */
            return( NULL );                   /* sorry.                       */
       }                                    
    else
        return( NULL );                       /* must be an 'X' or something.*/
   }

/********************************************/
/* make one entry into the permutation tree */
/********************************************/

void  tree_enter( TNODE *t, char *s, char *desc )
   {
    if ( verbose ) printf( "insert [%s]\n", s );
    if ( *s == '\0' )
       {
        t->desc = strdup( desc );
        /* printf( "inserted: %s\n", desc );  */
       }
    else
       {
        if ( t->ch[ind[*s]] == NULL )
           {
            /* printf( "need new node\n" );  */
            t->ch[ind[*s]] = new_tnode();
           }
        tree_enter( t->ch[ind[*s]], s+1, desc );
       }
   }

/***************************************************************************/
/* put a new pattern entry and its reverse-complement into the permutation */
/* tree                                                                    */
/***************************************************************************/

void  enter_pat( char *pat, int n_mis )
   {
    static  char  rpat[MAX_PAT_LEN];
    static  char  desc[MAX_STR_LEN];

    if ( verbose ) printf( "enter table [%s] %d\n", pat, n_mis );
    sprintf( desc, "%2d %s", n_mis, pat );
    tree_enter( t_root, pat, desc );
   }


                            /****************/
                            /* Permutations */
                            /****************/


/****************************************************/
/* convert all occurances of '.' to 'N' in string t */
/****************************************************/

void  unify_wildcards( char *t )
   {
    while ( *t != '\0' )
       {
        if ( *t == '.' ) *t = 'N';
        t++;
       }
   }


/* returns 1 if character c is a conserved nucleotide or ambiguity code. */
/* (uppercase indicates that its conserved                               */

int  is_conserved( int c )   
   {  return( isupper( c ) ); }


/*****************************************************************************/
/* permute() recursively generates all permuations of the pattern template,  */
/* entering each permution and its reverse complement into the permuation    */
/* tree.                                                                     */
/*                                                                           */
/*  char templ[] - pattern template string (command line argument), which    */
/*                 contains upper and/or lower case nucleotides and ambiguity*/
/*                 codes.                                                    */
/*  int max_mis  - this is the maximum number of total mismatches which are  */
/*                 are allowed.                                              */
/*  char s[]     - this is pre-allocated storage where the permuations are   */
/*                 "assembled"                                               */
/*  int  i       - recursion depth, or rather current position in templ[]    */
/*                 being considered.                                         */
/*  int  n_mis   - current number of accumulated mismatches at this          */
/*                 position, up to but not including s[i]                    */
/*                                                                           */
/*****************************************************************************/


void  permute( char templ[], int max_mis, char s[], int i, int n_mis )
   {
    char *matches;
    int   k, n;

    /* printf( "permute i = %d n_mis = %d\n", i, n_mis ); */

    if ( templ[i] == '\0' )   /* termination test - are we at end of templ?  */
       {                      /* if so, then s[] contains a permutation, so  */
        s[i] = '\0';          /* cap it off and then enter it and its rc into*/
        enter_pat( s, n_mis ); /* the permuation tree */
        if ( verbose ) printf( "%s %d\n", s, n_mis );
#ifdef SHOW_PERM_STATS
        num_perms++;
#endif
       }
    else /* if we have exhausted our mismatches at this point, or if */
         /* or if the template at this position is fixed (capitalized) */
        if ( n_mis >= max_mis  || is_conserved( templ[i] ) )  
           {    
            matches = allowed_matches[templ[i]];
            n = strlen( matches );
            for ( k = 0; k < n; k++ )   /* for each of the allowed possible */
               {                        /* nucleotides at this position */
                s[i] = matches[k];      /* insert it into s and then recurse */
                permute( templ, max_mis, s, i+1, n_mis );  /* do not update */
                                                           /* mismatch count*/
               }
           }
        else    /* otherwise, we're allowed at least one more mismatch at */
           {    /* this  position, so we'll generate four subpermutations */
            for ( k = 0; k < 4; k++ )
               {                                /* for each of A, C, G and T */
                s[i] = "ACGT"[k];               /* update s at this position */
                permute( templ, max_mis, s, i+1,  /* and then recurse, but   */
                         ( match("ACGT"[k],templ[i]) ? n_mis : (n_mis + 1) ) );
                                                  /* increment mismatch count*/
                                                  /* only if we mismatch   */
               }
           }
   }

void  generate_permutations( void )
   {
    static char  res[MAX_PAT_LEN];

    unify_wildcards( templ );

    init_tree();
    fill_allowed_matches(); 

    permute( templ, mis_thresh, res, 0, 0 ); 
#ifdef SHOW_PERM_STATS
    printf( "%d permutations\n", num_perms );
#endif
   }

/*********************************************************************/
/* open_file() opens a file or returns stdin if name is "-", or does */
/* error exit if file can't be opened                                */
/*********************************************************************/

FILE *open_file( char *name )
   {
    FILE *f;

    if ( strcmp( name, "-" ) == 0 )
        return( stdin );
    else
        if ( (f = fopen( name, "r" ) ) )
            return( f );
        else
           {
            perror( name );
            exit( errno );
           }
   }

/* close_file() closes the file unless its stdin */

void  close_file( FILE *f )
   {
    if ( f != stdin )
        fclose( f );
   }

                       /**********************/
                       /* comparision buffer */
                       /**********************/

static char  f_buff[MAX_BUFF_LEN+1];  /* top strand - possibly bisulf. mod. */
static char  r_buff[MAX_BUFF_LEN+1];  /* bottom strand   "        "     "   */
static char  u_buff[MAX_BUFF_LEN+1];  /* complement of top strand (after mod)*/
static char  v_buff[MAX_BUFF_LEN+1];  /* complement of bottom strand "    " */


void  init_buff( void )
   {
    size_t i;

    buff_len = neighbor_len + templ_len + neighbor_len;
    if ( buff_len > MAX_BUFF_LEN )
       {
        fprintf( stderr, "prosearch: buffer length %d exceeds max %d\n",
                          buff_len, MAX_BUFF_LEN );
        exit( 1 );
       }
    
    for ( i = 0; i < buff_len; i++ )
        r_buff[i] = f_buff[i] = 'X';
    r_buff[buff_len] = f_buff[buff_len] = '\0';

    if ( bisulfite_level > 0 )
       {
        for ( i = 0; i < buff_len; i++ )
            v_buff[i] = u_buff[i] = 'X';
        v_buff[buff_len] = u_buff[buff_len] = '\0';
       }

    if ( verbose )
        printf( "buffer initialized: %d + %d + %d = %d\n",
                     neighbor_len, templ_len, neighbor_len, buff_len );
   }


/* bisulfite_mod[c] contains the complementary nucleotide character for c */

char bisulfite_mod[] = {
                  /*          0     1     2     3     4     5     6     7   */
                  /* 000 */  'X',  'X',  'X',  'X',  'X',  'X',  'X',  'X',
                  /* 010 */  'X',  'X',  'X',  'X',  'X',  'X',  'X',  'X',
                  /* 020 */  'X',  'X',  'X',  'X',  'X',  'X',  'X',  'X',
                  /* 030 */  'X',  'X',  'X',  'X',  'X',  'X',  'X',  'X',
                  /* 040 */  'X',  'X',  'X',  'X',  'X',  'X',  'X',  'X',
                  /* 050 */  'X',  'X',  'X',  'X',  'X',  'X',  'X',  'X',
                  /* 060 */  'X',  'X',  'X',  'X',  'X',  'X',  'X',  'X',
                  /* 070 */  'X',  'X',  'X',  'X',  'X',  'X',  'X',  'X',
                  /* 100 */  'X',  'A',  'B',  'T',  'D',  'X',  'X',  'G',
                  /* 110 */  'H',  'X',  'X',  'K',  'X',  'M',  'N',  'X',
                  /* 120 */  'X',  'X',  'R',  'S',  'T',  'X',  'V',  'W',
                  /* 130 */  'X',  'Y',  'X',  'X',  'X',  'X',  'X',  'X',
                  /* 140 */  'X',  'a',  'b',  't',  'd',  'X',  'X',  'g',
                  /* 150 */  'h',  'X',  'X',  'k',  'X',  'm',  'n',  'X',
                  /* 160 */  'X',  'X',  'r',  's',  't',  'X',  'v',  'w',
                  /* 170 */  'X',  'y',  'X',  'X',  'X',  'X',  'X',  'X',
                 };

void  enter_base( int c )
   {
    size_t i;
    char   d;

    d = comp_char[c];

    if ( bisulfite_level > 0 )
       {
        c = bisulfite_mod[c];
        d = bisulfite_mod[d];

        memmove( u_buff+1, u_buff, buff_len-1 );
        u_buff[0] = comp_char[c];

        memmove( v_buff, v_buff+1, buff_len-1 );
        v_buff[buff_len-1] = comp_char[d];
       }
    memmove( f_buff, f_buff+1, buff_len-1 );   /* left shift f tag buffer */
    f_buff[buff_len-1] = c;

    memmove( r_buff+1, r_buff, buff_len-1 );
    r_buff[0] = d;
#ifdef DEBUG
    printf( "f %s\n", f_buff );
    printf( "u %s\n", u_buff );
    printf( "v %s\n", v_buff );
    printf( "r %s\n", r_buff );
    printf( "\n" );
#endif
   }


void  neighbor_output( char *desc, int pos, char *buf, char dir, 
                       char *filename )
   {
    static char pbuff[MAX_BUFF_LEN+1];
    static char left[MAX_NEIGHBOR_LEN+1];
    static char mat[MAX_PAT_LEN+1];
    static char fmt[MAX_STR_LEN+1];
    int   n_mis;

    /* TO FIX- don't scan dir from desc - get dir from dir,  */
    /*         use buf instead of buff, don't rc() */

    /* printf( "neigh %d [%s] [%s]\n", pos, buff, desc );  return; */

    sscanf( desc, "%d", &n_mis );
    strncpy( pbuff, buf, buff_len );
     
    strncpy( left, pbuff, neighbor_len );
    strncpy( mat,  pbuff + neighbor_len, templ_len );
    sprintf( fmt, "%%10d  %%c  %%s %%-%d.%ds %%-10s %%2d", 
                             templ_len, templ_len );
    printf( fmt,
             pos, dir, left, mat, pbuff+neighbor_len+templ_len, n_mis );
   }


void  lookup( char *buf, char *fbuf, char dir, int pos, 
              char *filename, char *hdr )
   {
    char        *rec;
    static char  pbuff[MAX_PAT_LEN+1];

    if ( (rec = tree_rlookup( t_root, buf + neighbor_len, templ_len,
                            amb_thresh ) ) )
       {
        if ( neighbor_len > 0 )
            neighbor_output( rec, pos, buf, dir,
                            (print_filenames ? filename : hdr ));
        else
           {
            strncpy( pbuff, fbuf+neighbor_len, templ_len );
            printf( "%s %c %s %10d", pbuff, dir, rec, pos );
           }
        if ( print_headers )
            printf( " %s\n", hdr );    /* hdr has a \n in it */
        else
            putchar( '\n' );
       }
   }

/*****************************************************/
/* scan_file() - open file and process its sequences */
/*****************************************************/

void  scan_file( char *filename )
   {
    FILE        *f;
    int          c;
    int          last;
    static char  hdr[MAX_HDR_LEN+1];
    int          l;     /* header length */
    int          pos;   /* position within string */

    if ( verbose )
        printf( "file %s\n", filename );
    f = open_file( filename );
    last = '\n';
    if ( verbose ) 
        printf( "templ_len is %d\n", templ_len );
    while ( (c = fgetc( f ) ) != EOF )
       {
        if ( last == '\n' && c == '>' )
           {
            fgets( hdr, MAX_HDR_LEN, f );
            if ( verbose )
                printf( "seq: %s", hdr );
            l = strlen( hdr );
            if ( l > 0 && hdr[l-1] == '\n' )       /* behavior of fgets() */
                hdr[l-1] = '\0';
            init_buff();
            pos = 1 - (templ_len + neighbor_len);  /* counting from one */
            last = '\n';
           }
        else if ( isalpha( c ) )
           {
            enter_base( c );
            ++pos;
            /* if ( rec = tree_lookup( buff + neighbor_len, templ_len ) ) */
            lookup( f_buff, f_buff, 'f', pos, filename, hdr );
            lookup( r_buff, f_buff, 'r', pos, filename, hdr );
            if ( bisulfite_level > 0 )
               {
                lookup( u_buff, f_buff, 'u', pos, filename, hdr );
                lookup( v_buff, f_buff, 'v', pos, filename, hdr );
               }
            last = c;
           }
        else
            last = c;
       }
    close_file( f );
   }



                              /****************/
                              /* Main Program */
                              /****************/


int  main( int argc, char **argv )
   {
    int     nfiles;
    char  **filenames;
    int     i;

    parse_args( argc, argv, &nfiles, &filenames );
    templ_len = strlen( templ );

    if ( verbose )
        printf( "mismatch threshold: %d\n", mis_thresh );

    generate_permutations();

    for ( i = 0; i < nfiles; i++ )
        scan_file( filenames[i] );
   }
