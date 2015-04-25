/* Program:       gsts2                                                      */
/* Programmer:    Sean R. McCorkle                                           */
/*                Biology Department, Brookhaven National Laboratory         */
/* Language:      C                                                          */
/*                                                                           */
/* Description:   C version of original perl version of gsts                 */
/*                                                                           */
/* Usage:         gsts [-chvFIV] [-a<seq>] [-e<n>] [-f<seq>] [seq-file ... ] */
/*                                                                           */
/*                where [seq-files] are in FASTA format.  The name "-" means */
/*                stdin.  stdin is scanned it no filemames are specified.    */
/*                                                                           */
/* Options:       -a<seq> use <seq> as anchor enzyme sequence (def. CATG)    */
/*                -c      treat all sequences as circular                    */
/*                -e<n>   tag extent <n>, (default 17)                       */
/*                -f<seq> use <seq> as fragmenting enz. seq. (def. ACTAGT)   */
/*                -h      print help message                                 */
/*                -l<seq> append linker sequence <seq> to short tags         */
/*                        (default GGATCCGAAGGGGTTCG)                        */
/*                -n<n>   number sequences starting from <n>                 */
/*                -v      verbose mode                                       */
/*                -F      print out fragment information, rather than tag    */
/*                -I      print out internal tags as well as true GSTs       */
/*                -S      print out anchor site to 3' fragment end sequences */
/*                        (fasta format).  Doesn't include linker            */
/*                -E      print excluded sequence from each fragment         */
/*                        (between -S sequences)
/*                -T      print fragment & anchor totals                     */
/*                -V      print version                                      */
/*                                                                           */
/* Caveats:       Sequence characters other than standard nucleotide and     */
/*                ambiguity codes will not be detected and will cause        */
/*                undefined behavior.  GIGO                                  */
/*                                                                           */
/*                Beware of -I (internals) behavior when there are no        */
/*                fragmenting enzymes in the sequence at all - default is    */
/*                to print nothing - probably this should be changed         */
/*                                                                           */
/* Compiling:     cc -O -o gsts2 gsts2.c                                     */
/*                                                                           */
/*                  - usually does the job on most unix platforms            */
/*                    (including MacOSX)                                     */
/*                                                                           */
/* NOTE TO SELF: absolute 5' positions are off for non-circular fragment 0   */
/*               also: rationalize the position scheme once and for all      */
/*                                                                           */
/* Notes:         What is lacking, but eventually planned:                   */
/*                    1) restriction enzyme sequences can't contain          */
/*                       ambiguity codes (or be regexps).                    */
/*                    2) ambiguity codes in sequences do not match properly  */
/*                    3) short tags are properly handled, but not internal   */
/*                       tags with short lengths between anchor sites        */
/*                                                                           */
/* $Id: gsts2.c,v 0.7 2003/03/03 15:43:05 mccorkle Exp mccorkle $            */
/*                                                                           */
/*****************************************************************************/

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifndef SOLARIS        /* Sun's are really a pain in the neck sometimes */
#include <unistd.h>    /* for getopt(), BSD, Irix, Linux */
#endif

       /* various globals set by command line args and options */

char *anch_enz = "CATG";   /* anchor enzyme seq., set by -a, default NlaIII */
char *frag_enz = "ACTAGT"; /* fragmnt enzyme seq., set by -f, default SpeI  */
int   tag_ext  = 17;       /* tag extent, set by -e, default for MmeI       */
int   totals   = 0;        /* show fragment & anchor totals, set by -T      */
int   verbose  = 0;        /* verbose mode, set by -v                       */
int   internals= 0;        /* show internal tags, set by -I                 */
int   fragments= 0;        /* show fragment info, set by -F                 */
int   extensions=0;        /* show anchor->fragment extensions, set by -S   */
char *linker   = "GGATCCGAAGGGGTTCG";         /* linker sequence, set by -l */


#define  MAX_TAG_EXT     100
#define  FASTA_LINE_LEN   60

int  anch_enz_len;
int  frag_enz_len;

       /* state variables */

int  pos;                  /* current position in genome sequence.  set and  */
                           /* incremented by init_cbuff() and insert()       */

int  seq_num;              /* current sequence number.  initialized in parse_*/
                           /* args() (-n), and incremented in new_sequence() */
                           

int  frag_num;             /* current fragment number;  0 indicates region   */
                           /* between 5' end of genome sequence and first    */
                           /* fragmenting enzyme site.  (re)set by           */
                           /* handle_last_frag() and incremented in          */
                           /* handle_fragment().                             */

int  last_frag_pos;        /* position, in genome sequence, just after the 3'*/
                           /* end of the last frag. enzyme site encountered. */
                           /* (re)set in handle_last_frag() and incremented  */
                           /* in handle_fragment()                           */

/* we have a tag list for internal tags in a fragment.  */

#define  MAX_TAGS           20000
char    *ftags[MAX_TAGS];   /* forward tags */
char    *rtags[MAX_TAGS];   /* backward tags */
int      tag_pos[MAX_TAGS]; /* position of tag (left end of anchor site) */
int      ntags;

/* duplicate for fragment 0 (5' initial fragment of each sequence */

char    *ftags_zero[MAX_TAGS];   /* forward tags */
char    *rtags_zero[MAX_TAGS];   /* backward tags */
int      tag_pos_zero[MAX_TAGS]; /* ... etc */
int      ntags_zero;
int      frag_len_zero;

/* character string storage for 3' extensions */

#define  MAX_EXT_SEQ   2000000
char*    ext_next;                    /* global pointer */
char     ext_buff[MAX_EXT_SEQ+1];     /* these are null terminated */
char     r_extension[MAX_EXT_SEQ+1];
int      r_ext_set;                   /* boolean flag, cleared by reset_     */
                                      /* extensions(), set by handle_anchor()*/
char*    f_ext_zero;                  /* what a kludge.  fix this. */     

/* counters */

int  num_anchs = 0;         /* counts total number of anchor sites */
int  num_frags = 0;         /* counts total number of fragment sites */
int  num_frags_w_gsts = 0;  /* total number of fragments with GSTs */
int  num_frags_wo_gsts = 0; /* total number of fragments without GSTs */

static char const rcsid[] =
    "$Id: gsts2.c,v 0.7 2003/03/03 15:43:05 mccorkle Exp mccorkle $";
static char rcsvers[] = "$Revision: 0.7 $";


/* version() - print program name and version */

void  version( void )
   {
    char *v;
    for ( v = rcsvers; *v && *v != ' '; v++ ) 
       ;
    v++;
    v[strlen(v)-1] = '\0';
    printf( "\ngsts2        version %s\n", v );
   }



/* version() - print program name, version and help */

void  help( void )
   {
    version();
    printf( " \n\
Usage:       gsts [-chvFIV] [-a<seq>] [-e<n>] [-f<seq>] [seq-file ... ]   \n\
                                                                          \n\
             where [seq-files] are in FASTA format.  The name \"-\" means \n\
             stdin.  stdin is scanned it no filemames are specified.      \n\
                                                                          \n\
Options:     -a<seq> use <seq> as anchor enzyme sequence (def. CATG)      \n\
             -c      treat all sequences as circular                      \n\
             -e<n>   tag extent <n>, (default 17)                         \n\
             -f<seq> use <seq> as fragmenting enz. seq. (def. ACTAGT)     \n\
             -h      print help message                                   \n\
             -l<seq> append linker sequence <seq> to short tags           \n\
                     (default GGATCCGAAGGGGTTCG)                          \n\
             -n<n>   number sequences starting from <n>                   \n\
             -v      verbose mode                                         \n\
             -F      print out fragment information, rather than tag      \n\
             -I      print out internal tags as well as true GSTs         \n\
             -S      print out anchor site to 3' fragment end sequences   \n\
                     (fasta format).  Doesn't include linker              \n\
             -T      print fragment & anchor totals                       \n\
             -V      print version                                        \n\
                                                                          \n\
" );

   }

void  notyet( char opt ) 
   { fprintf( stderr, "sorry: %c option is not implemented yet\n", opt ); 
     exit(1);
   }

/* uc( str ) - convert null-terminated str into uppercase, return ptr to s */

char *uc( char *s )
   {
    char *t;
    
    for ( t = s; *t; t++ )
        *t = toupper( *t );
    return( s );
   }

/**************************************************************************/
/* parse_ags() - read command line arguments and options and set globals. */
/* returns remaining arguments as nfiles and filenames array              */
/* note to self: the invocations of atoi() and strdup() could probably be */
/* made safer with more checks                                            */
/**************************************************************************/

/* note to self: put in length checks on -a, -f, -l, etc */

void  parse_args( int argc, char **argv, int *nfiles, char ***filenames )
   {
    extern char *optarg;
    extern int   optind;
    int          c;
    static char *def_files[] = { "-", "" };

    seq_num = 1;
    while ( (c = getopt( argc, argv, "a:ce:f:hl:n:vFISTV" ) ) != -1 )
        switch( c )
           {
            case  'a':   anch_enz = uc( strdup( optarg ) );
                         break;
            case  'c':   notyet( 'c' );
                         break;
            case  'e':   tag_ext = atoi( optarg );
                         if ( tag_ext > MAX_TAG_EXT || tag_ext < 0 )
                            {
                             fprintf( stderr, 
                               "tag extent must be non-negative & < %d\n",
                               MAX_TAG_EXT );
                             exit( 1 );
                            }
                         break;
            case  'f':   frag_enz = uc( strdup( optarg ) );
                         break;
            case  'h':   help();
                         exit( 1 );
                         break;
            case  'l':   linker = uc( strdup( optarg ) );
                         break;
            case  'n':   seq_num = atoi( optarg );
                         break;
            case  'v':   verbose = 1;
                         break;
            case  'F':   fragments = 1;
                         break;
            case  'I':   internals = 1;
                         break;
            case  'S':   extensions = 1;   /* future: handle conflict w/F,I..*/
                         break;
            case  'T':   totals = 1;
                         break;
            case  'V':   version();
                         exit( 0 );
            default:     help();
                         exit( 1 );
           }

    seq_num--;            /* for convenience */

    argc -= optind;
    argv += optind;
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



void  init_tag_list( tag_size )
   {
    int i;

    if ( verbose ) printf( "initializing tag queue...\n" );
    for ( i = 0; i < MAX_TAGS; i++ )
       {
        if ( ! (ftags[i] = calloc( tag_size+1, 1 ))  ||
             ! (rtags[i] = calloc( tag_size+1, 1 ))  ||
             ! (ftags_zero[i] = calloc( tag_size+1, 1 ))  ||
             ! (rtags_zero[i] = calloc( tag_size+1, 1 ))
           )
           {
            fprintf( stderr, 
             "Fatal error: unable to allocate %d bytes for tag queue pos %d\n",
                 tag_size, i );
            exit( 1 );
           }
       }
    if ( verbose ) printf( "...done\n" );
    ntags = 0;
    ntags_zero = 0;
   }

void  insert_tags( char *ftag, char *rtag, int tag_len )
   {
    if ( pos - last_frag_pos < 0 )  /* don't insert if anchor site */
        return;                     /* overlaps with fragment site */
    tag_pos[ntags] = pos - last_frag_pos;
    strncpy( ftags[ntags], ftag, tag_len );
    strncpy( rtags[ntags++], rtag, tag_len );
    if ( ntags >= MAX_TAGS )
       {
        fprintf( stderr, 
         "Fatal error: number of internal tags has exceeded maximum of MAX_TAGS = %d\n",
                 MAX_TAGS );
        fprintf( stderr,"try increasing value of MAX_TAGS and recompiling\n");
        exit( 1 );
       }
   }

void  copy_tags_zero( void )
   {
    int  i;

    for ( i = 0; i < ntags; i++ )
       {
        tag_pos_zero[i] = tag_pos[i];
        strncpy ( ftags_zero[i], ftags[i], tag_ext );
        strncpy ( rtags_zero[i], rtags[i], tag_ext );
       }
    ntags_zero = ntags;
    frag_len_zero = pos;
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
            exit( 1 );
           }
   }

/* close_file() closes the file unless its stdin */

void  close_file( FILE *f )
   {
    if ( f != stdin )
        fclose( f );
   }


/* circular buffer is just a shifting q for now.  Improvements later */

#define MAX_CBUFF  100 
char   cbuff[MAX_CBUFF];

int    enz_pos;
int    cbuff_len;

void  init_cbuff( void )
   {
    if ( verbose ) printf( "initializing cbuff..\n" );
    cbuff_len = 2 * tag_ext + anch_enz_len;  /* could check against MAX_CBUFF*/
    memset( cbuff, 'X', cbuff_len );
    cbuff[cbuff_len] = '\0';
    pos = -( anch_enz_len + tag_ext );
   }


void  insert( char c )
   {
    if ( verbose ) printf( "insert %c ", c );
    memmove( cbuff, cbuff+1, cbuff_len-1 );
    cbuff[cbuff_len - 1] = toupper( c );
    pos++;
    if ( verbose ) printf( "%10d %s\n", pos, cbuff );
   }

int  is_anchor_site( void )
   { return( strncmp( cbuff + tag_ext, anch_enz, anch_enz_len ) == 0 ); }
    

int  is_fragmenting_site( void )
   { return( strncmp( cbuff + tag_ext, frag_enz, frag_enz_len ) == 0 ); }

/*returns the current character (i.e. character in sequence at position pos) */

char  curr_char( void )
   { return( cbuff[tag_ext] ); }

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


/* rc( s, n, r ) - put reverse complment of s (length n) into r */

void  rc( char *s, int n, char *r )
   {
    while ( n > 0 )
        *r++ = comp_char[s[--n]];
    *r = '\0';
   }

/* write out the null-terminated sequence in seq as a FASTA body (no header) */

void  output_seq( char *seq )
   {
    char *s;
    int   i = 0;

    for ( s = seq;  *s != '\0';  s++ )
       {
        putchar( *s );
        if ( ( ++i % FASTA_LINE_LEN ) == 0 )
            putchar( '\n' );
       }
    if ( i % FASTA_LINE_LEN ) putchar( '\n' );
   }

void  reset_ext_buff( void )
   {
    if ( verbose ) printf( "Reseting extension buff\n" );
    ext_next = ext_buff;
   }

void  reset_extensions( void )
   {
    if ( verbose ) printf( "Reseting extensions\n" );
    r_ext_set = 0;
    reset_ext_buff();
   }

void  add_ext_char( char c )
   {
    if ( ext_next >= ext_buff + MAX_EXT_SEQ )
       {
        fprintf(stderr,"extension sequence too long (> %d)\n", MAX_EXT_SEQ );
        fprintf(stderr,"  (try different enzymes or increase MAX_EXT_SEQ)\n" );
        exit( 1 );
       }
    *ext_next++ = c;
   }

void  print_f_ext( void )
   {
    if ( ntags > 0 )
        {
         add_ext_char( '\0' );
         output_seq( ext_buff + anch_enz_len - 1 );
        }
   }

void  print_r_ext( void )
   {
    if ( ntags > 0 )
        output_seq( r_extension );
   }

void  save_f_ext_zero( void )
   {
    if ( ntags > 0 )
        {
         add_ext_char( '\0' );
         f_ext_zero = strdup( ext_buff + anch_enz_len - 1 );
        }
   }


void  handle_anchor( void )
   {
    static char rtag[MAX_TAG_EXT];
    int    n;

    if ( verbose ) printf( "     handling anchor\n" );

    rc( cbuff, tag_ext, rtag );
    insert_tags( cbuff + tag_ext + anch_enz_len, rtag, tag_ext );
    if ( extensions )
       {
        if ( ! r_ext_set )
           {
            n = (ext_next - ext_buff) - (frag_enz_len - 1);
            if ( n < 0 ) 
                n = 0;
            rc( ext_buff + frag_enz_len - 1, n,  r_extension );
            r_ext_set = 1;
           }
        reset_ext_buff();
       }
    num_anchs++;   /* incr. total count seen */
   }

void  print_tags( char *ftags[], char *rtags[], int tag_pos[], int ntags,
                  int  frag_num, int  frag_len )
   {
    int          i;
    static char *flags="";
    int          j, k;
   
    for ( i = 0; i < ntags; i++ )
       {
        if ( frag_len - (tag_pos[i] + anch_enz_len) < tag_ext ) 
           {
            flags = "S";
            for ( j = frag_len - (tag_pos[i] + anch_enz_len), k = 0;
                  j < tag_ext;
                  j++, k++ )
                (ftags[i])[j] = linker[k];
           }
        else
            flags = "";
        printf( "%-24.24s %3d f %5d %8d %4d %8d %4d  %10d %s\n",
                  ftags[i], seq_num, frag_num, 
                 (frag_len - tag_pos[i]), ntags - i, tag_pos[i], i+1, 
                 tag_pos[i] + last_frag_pos, flags );

        if ( tag_pos[i] < tag_ext ) 
           {
            flags = "S";
            for ( j = tag_pos[i], k = 0;  j < tag_ext;   j++, k++ )
                (rtags[i])[j] = linker[k];
           }
        else
            flags = "";

        printf( "%-24.24s %3d r %5d %8d %4d %8d %4d  %10d %s\n",
                 rtags[i], seq_num, frag_num, 
                 (tag_pos[i] + anch_enz_len), i+1, 
                 (frag_len - (tag_pos[i] + anch_enz_len)), ntags - i, 
                 tag_pos[i] + last_frag_pos, flags  );
        flags = "";
       }
   }

void  print_ftag( char *ftags[], int tag_pos[], int ntags,
                  int  frag_num, int  frag_len )
   {
    int          i;
    static char *flags="";
    int          j, k;

    if ( ntags > 0 )
       {
        i = ntags - 1;
        if ( frag_len - (tag_pos[i] + anch_enz_len) < tag_ext ) 
           {
            flags = "S";
            for ( j = frag_len - (tag_pos[i] + anch_enz_len), k = 0;
                  j < tag_ext;
                  j++, k++ )
                 (ftags[i])[j] = linker[k];
            }
        else 
            flags = "";
        if ( extensions ) printf( ">" );
        printf( "%-24.24s %3d f %5d %8d %4d %8d %4d  %10d %s\n",
                  ftags[i], seq_num, frag_num, 
                 (frag_len - tag_pos[i]), ntags - i, tag_pos[i], i+1, 
                 tag_pos[i] + last_frag_pos, flags );
        }
   }

void  print_rtag( char *rtags[], int tag_pos[], int ntags,
                  int  frag_num, int  frag_len )
   {
    static char *flags="";
    int          j, k;
   
    if ( ntags > 0 )
       {
        if ( tag_pos[0] < tag_ext ) 
           {
            flags = "S";
            for ( j = tag_pos[0], k = 0;  j < tag_ext;   j++, k++ )
                (rtags[0])[j] = linker[k];
           }
        else
            flags = "";

        if ( extensions ) printf( ">" );
        printf( "%-24.24s %3d r %5d %8d %4d %8d %4d  %10d %s\n",
                 rtags[0], seq_num, frag_num, 
                 (tag_pos[0] + anch_enz_len), 1, 
                 (frag_len - (tag_pos[0] + anch_enz_len)), ntags, 
                 tag_pos[0] + last_frag_pos, flags  );
        flags = "";
       }
   }


void  print_frags( int frag_num, int frag_len )
   {
    printf( "%8d %8d %8d %8d %8d\n", 
            seq_num, frag_num, last_frag_pos, frag_len, ntags );
   }

void  handle_fragment( void )
   {
    int frag_len = pos - last_frag_pos;  /* calculate fragment length */

    if ( verbose ) printf( "     handling fragment\n" );
    /* remove last tag if anchor enz. site overlaps with frag. enz. site */
    if ( (frag_len - (tag_pos[ntags-1] + anch_enz_len)) < 0 )
        ntags--;

    if ( frag_num > 0 )
       {
        if ( fragments )
            print_frags( frag_num, frag_len );
        else if ( internals )
            print_tags( ftags, rtags, tag_pos, ntags, frag_num, frag_len );
        else
           {
            print_ftag( ftags, tag_pos, ntags, frag_num, frag_len );
            if ( extensions )
                print_f_ext();
            print_rtag( rtags, tag_pos, ntags, frag_num, frag_len );
            if ( extensions )
                print_r_ext();
           }
        if ( ntags > 0 )
            num_frags_w_gsts++;
        else
            num_frags_wo_gsts++;
       }
    else
       {
        copy_tags_zero();  /* note: this is handling overlapping anchor & */
                           /* and fragment sites in frag 0                */
        if ( extensions )
            save_f_ext_zero();
       }
    if ( extensions )
        reset_extensions();
    frag_num++;
    last_frag_pos = pos + frag_enz_len;
    ntags = 0;                          /* resets tag list */
    num_frags++;                        /* incr total count */
   }


/* note to self:  this has a lot of kind of redundancy - maybe there's */
/* a better way?                                                       */

void  handle_last_frag( void )
   {

    if ( verbose ) printf( "     handling last fragment\n" );


    if ( frag_num > 0 )   /* default is to not print if no fragment enzyme   */
       {                  /* sites are found - maybe this could be an option?*/
        if ( fragments )
            print_frags( 0, frag_len_zero );
        else if ( internals )
           {
            /* note: frag 0 and last frag overlapping anchor enz. and frag. */
            /* enz. sites are already taken care of                         */
            print_tags( ftags_zero, rtags_zero, tag_pos_zero, ntags_zero, 0,
                        frag_len_zero );
            print_tags( ftags, rtags, tag_pos, ntags, frag_num,
                         pos - last_frag_pos );
           }
        else
           {
            print_ftag(ftags_zero, tag_pos_zero, ntags_zero, 0,frag_len_zero);
            if ( extensions && ntags_zero > 0 )
                output_seq( f_ext_zero );
            print_rtag( rtags, tag_pos, ntags, frag_num, pos - last_frag_pos );
            if ( extensions )
                print_r_ext();
           }
        if ( ntags > 0 )
            num_frags_w_gsts++;
        else
            num_frags_wo_gsts++;
       }
    if ( extensions )
        reset_extensions();
    frag_num = 0;
    last_frag_pos = 0;
    ntags = 0;                          /* resets tag list */
   }

/* handle( c ) - handle one input character  */

void  handle( char c )
   {
    if ( isalpha( c ) )
       {                              /* if its a nucleotide,    */
        insert( c );                  /* then put it into cbuff  */
        if ( is_anchor_site() )       /* and check to see if its */
            handle_anchor();          /* one of our enzyme seqs  */
        else if ( is_fragmenting_site() ) /* assume that anchor &*/
            handle_fragment();        /* frag sites are mutually */
                                      /* exclusive               */
        else if ( extensions )
            add_ext_char( curr_char() );
       }
   }

/* this runs out the last part of the genome sequence through the cbuff */
/* to the matching position.                                            */

void  flush_last_seq( void )
   {
    int   i;

    for ( i = 0; i < tag_ext + anch_enz_len; i++ )
       {
        if ( verbose ) printf( "flushing %d\n", i );
        handle( 'X' );
       }
   }


void  new_sequence( FILE *f )
   {
    int  c;

    if ( verbose ) printf( "new sequence\n" );
    flush_last_seq(); /* flush out last bit of buffer */
    handle_last_frag();
    init_cbuff();     /* new sequence => starting with fresh cbuff and fresh */
                      /* state variables */
    while ( (c = fgetc( f )) != EOF && c != '\n' )
        if ( verbose ) putchar( c );
    if ( verbose ) printf( "  char is now [%d]\n", c );
    f_ext_zero = "";
    seq_num++;
   }


void  show_totals( void )
   {
    printf( "total fragment sites: %d\n", num_frags );
    printf( "total anchor sites: %d\n", num_anchs );
    printf( "total fragments containing GSTs: %d\n", num_frags_w_gsts );
    printf( "total fragments with no GSTs: %d\n", num_frags_wo_gsts );

   }


                                 /****************/
                                 /* Main program */
                                 /****************/

main ( int argc, char **argv )
   {
    int     nfiles;
    char  **filenames;
    int     i;
    FILE   *f;
    int     c;
    char    last_char;


    parse_args( argc, argv, &nfiles, &filenames );

    anch_enz_len = strlen( anch_enz );
    frag_enz_len = strlen( frag_enz );

    init_tag_list( tag_ext );
    init_cbuff();
    if ( extensions )
        reset_extensions();
    f_ext_zero = "";

    if ( verbose )
       {
        printf( "anchor:     %s\n", anch_enz );
        printf( "anchor len: %d\n", anch_enz_len );
        printf( "fragmenter: %s\n", frag_enz );
        printf( "frag len:   %d\n", frag_enz_len );
        printf( "linker:     %s\n", linker );
        printf( "tag extent: %d\n", tag_ext );
        printf( "nfiles:     %d\n", nfiles );
        for ( i = 0; i < nfiles; i++ )
           printf( "               %s\n", filenames[i] );
       }
    
    for ( i = 0; i < nfiles; i++ )
       {
        if ( verbose )
            printf( "File:   %s\n", filenames[i] );
        f = open_file( filenames[i] );
        last_char = '\n';
        while ( (c = fgetc( f )) != EOF )
           {
            if ( c == '>' && last_char == '\n' )
               {
                new_sequence( f );
                last_char = '\n';
               }
            else 
               {
                handle( c );
                last_char = c;
               }
           }
        close_file( f );
       }
    flush_last_seq();
    handle_last_frag();
    if ( totals )
        show_totals();
   }
