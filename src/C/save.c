/* Program:     prosearch.c                                                  */
/* Programmer:  Sean R. McCorkle                                             */
/*              Biology Dept., Brookhaven National Laboratory                */
/* Language:    C                                                            */
/*                                                                           */
/* Description: Search DNA sequence files (FASTA format) for candidate       */
/*              promoter sequences, with adjustable mismatch threshold       */
/*                                                                           */
/* Usage:       prosearch [-m<n>] [-hvV] <pat> [<file> [...]]                */
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
/* Options:     -m<n>  accept up to <n> mismatches in LOWER case bases in    */
/*                     <pat>.  Default is 0 (exact matches).                 */
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
/*               Ambiguity codes in the input sequence files are NOT handled */
/*               properly.  They are ignored (for now).                      */
/*                                                                           */
/* $Id: prosearch.c,v 0.6 2003/11/10 20:52:16 mccorkle Exp $        */
/*                                                                           */
/*****************************************************************************/

#include <ctype.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

static char proc_rcs_id[] =
    "$Id: prosearch.c,v 0.6 2003/11/10 20:52:16 mccorkle Exp $";

#define MAX_HDR_LEN 256
#define MAX_PAT_LEN 256
#define MAX_STR_LEN 256

#define NULL        ((void *) 0)

char  bases[] = "ACGT";
/* char  templ[] = "ttcaGCacc.cGGacagc.cc"; */
char  templ[MAX_PAT_LEN+1];
int   templ_len;


int   verbose = 0;      /* set by -v */
int   mis_thresh = 0;   /* set by -m */



void  version( void )
   {
    char *v;
    v = strdup("$Revision: 0.6 $"+11 );
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
Usage:       prosearch [-m<n>] [-hvV] <pat> [<file> [...]]                  \n\
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
Options:     -m<n>  accept up to <n> mismatches in LOWER case bases in      \n\
                    <pat>.  Default is 0 (exact matches).                   \n\
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


void  parse_args( int argc, char **argv, int *nfiles, char ***filenames )
   {
    extern char *optarg;
    extern int   optind;
    int          c;
    static char *def_files[] = { "-", "" };

    while ( (c = getopt( argc, argv, "hm:vV" ) ) != -1 )
        switch( c )
           {
            case  'h':   help();
                         exit( 0 );
            case  'm':   mis_thresh = atoi( optarg );
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

char *tree_lookup( char *s )
   {
    TNODE *t;

    if ( verbose )
        printf( "lookup [%s]\n", s );
    t = t_root;
    while ( *s != '\0' && ind[*s] >= 0 && t->ch[ind[*s]] !=NULL )
        t = t->ch[ind[*s++]];

    if ( *s == '\0' )
        return( t->desc );
    else if ( ind[*s] < 0 )
        return( NULL );       /* 'X' or anything else in s fails */
    else
        return( NULL );
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
    sprintf( desc, "f %2d %s", n_mis, pat );
    tree_enter( t_root, pat, desc );

    rc( rpat, pat );
    sprintf( desc, "r %2d %s", n_mis, pat );
    tree_enter( t_root, rpat, desc );
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


/* allowed_matches maps each nucleotide and ambiguity code to an array  */
/* of all the nucleotides (ACGT) that will match it.                    */

char *allowed_matches[128]; 

/* this initializes the allowed_matches table */

void  fill_allowed_matches( void )
   {
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
    size_t       num_perms = 5000000;
    static char  res[MAX_PAT_LEN];

    unify_wildcards( templ );

    init_tree();
    fill_allowed_matches(); 

    permute( templ, mis_thresh, res, 0, 0 ); 
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

static char  buff[MAX_PAT_LEN+1];

void  init_buff( void )
   {
    size_t i;
    for ( i = 0; i < templ_len; i++ )
        buff[i] = 'X';
    buff[templ_len] = '\0';
   }

void  enter_base( int c )
   {
    size_t i;

    for ( i = 1; i < templ_len; i++ )   /* memmove would be better */
        buff[i-1] = buff[i];
    buff[templ_len-1] = c;
   }


/*****************************************************/
/* scan_file() - open file and process its sequences */
/*****************************************************/

void  scan_file( char *filename )
   {
    FILE        *f;
    int          c, last;
    static char  hdr[MAX_HDR_LEN+1];
    int          pos;   /* position within string */
    char        *rec;

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
            init_buff();
            pos = -templ_len;
            last = '\n';
           }
        else if ( isalpha( c ) )
           {
            enter_base( c );
            ++pos;
            if ( rec = tree_lookup( buff ) )
                printf( "%s %s %10d %s",  buff, rec, pos, hdr );
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


main( int argc, char **argv )
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
