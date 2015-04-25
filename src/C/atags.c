/* Program:       atags                                                       */
/* Programmer:    Sean R. McCorkle  (mccorkle@bnl.gov)                        */
/*                Biology Department, Brookhaven National Laboratory          */
/* Language:      C                                                           */
/*                                                                            */
/* Description:   Read genome sequence, print all GST tags derived from a     */
/*                a specified restriction enzyme                              */
/*                                                                            */
/* Usage:         atags [-hFSvV] [-a<seq>] [-e<n>] [seq-file ... ]            */
/*                                                                            */
/*                where [seq-files] are in FASTA format.  The name "-" means  */
/*                stdin.  stdin is scanned it no filemames are specified.     */
/*                                                                            */
/* Options:       -a<seq> use <seq> as anchor enzyme sequence (def. CATG)     */
/*                -e<n>   tag extent <n>, (default 17)                        */
/*                -h      print help message                                  */
/*                -l<seq> append linker sequence <seq> to short tags          */
/*                        (default GGATCCGAAGGGGTTCG)                         */
/*                -n<n>   number sequences starting at <n>                    */
/*                -F      Show fragment numbers                               */
/*                -S      Show sequence names instead of numbers              */
/*                -v      verbose mode                                        */
/*                -V      print version                                       */
/*                                                                            */
/* Caveats:       Sequence characters other than standard nucleotide and      */
/*                ambiguity codes will not be detected and will cause         */
/*                undefined behavior.  GIGO                                   */
/*                                                                            */
/* Compiling:     cc -O -o atags atags.c                                      */
/*                                                                            */
/*                  - usually does the job on most unix platforms             */
/*                    (including MacOSX)                                      */
/*                                                                            */
/* Notes:         What is lacking, but eventually planned:                    */
/*                    1) restriction enzyme sequences can't contain           */
/*                       ambiguity codes (or be regexps).                     */
/*                    2) ambiguity codes in sequences do not match properly   */
/*                    3) shorts are proprely handled, but not internal        */
/*                       tags with short lengths between anchor sites         */
/*                                                                            */
/* $Id: atags.c,v 1.1 2005/06/13 21:05:27 mccorkle Exp mccorkle $             */
/*                                                                            */
/******************************************************************************/

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifndef SOLARIS        /* Sun's are really a pain in the neck sometimes */
#include <unistd.h>    /* for getopt(), BSD, Irix, Linux */
#endif


       /* various globals set by command line args and options */

char *anch_enz = "CATG";   /* anchor enzyme seq., set by -a, default NlaIII */
int   tag_ext  = 17;       /* tag extent, set by -e, default for MmeI       */
int   verbose  = 0;        /* verbose mode, set by -v                       */
int   internals= 0;        /* show internal tags, set by -I                 */
char *linker   = "GGATCCGAAGGGGTTCG";         /* linker sequence, set by -l */
int   seqnames = 0;        /* set by -S */
int   fragnumbers = 0;     /* set by -F */

#define  MAX_TAG_EXT    10000
#define  MAX_HDR_LEN      256
#define  MAX_ANCH_LEN     100

char header[MAX_HDR_LEN+1];

int  anch_enz_len;

       /* state variables */

int  pos;                  /* current position in genome sequence.  set and  */
                           /* incremented by init_cbuff() and insert()       */

int  seq_num = 0;          /* current sequence number.  possibly set by -n,  */
                           /* and incremented in new_sequence()              */

int  frag_num = 0;



static char const rcsid[] =
    "$Id: atags.c,v 1.1 2005/06/13 21:05:27 mccorkle Exp mccorkle $";
static char rcsvers[] = "$Revision: 1.1 $";


/* version() - print program name and version */

void  version( void )
   {
    char *v;
    for ( v = rcsvers; *v && *v != ' '; v++ ) 
       ;
    v++;
    v[strlen(v)-1] = '\0';
    printf( "\natags        version %s\n", v );
   }



/* version() - print program name, version and help */

void  help( void )
   {
    version();
    printf( " \n\
Usage:       atags [-hFSvV] [-a<seq>] [-e<n>] [-l<seq>] [seq-file ... ]     \n\
                                                                          \n\
             where [seq-files] are in FASTA format.  The name \"-\" means \n\
             stdin.  stdin is scanned it no filemames are specified.      \n\
                                                                          \n\
Options:     -a<seq> use <seq> as anchor enzyme sequence (def. CATG)      \n\
             -e<n>   tag extent <n>, (default 17)                         \n\
             -h      print help message                                   \n\
             -l<seq> append linker sequence <seq> to short tags           \n\
                     (default GGATCCGAAGGGGTTCG)                          \n\
             -n<n>   number sequences starting at <n>                     \n\
             -F      Show fragment numbers                                \n\
             -S      Show sequence names instead of numbers               \n\
             -v      verbose mode                                         \n\
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

/* note to self: put in length checks on -a, -l, etc */

void  parse_args( int argc, char **argv, int *nfiles, char ***filenames )
   {
    extern char *optarg;
    extern int   optind;
    int          c;
    static char *def_files[] = { "-", "" };

    while ( (c = getopt( argc, argv, "a:e:Fhl:n:SvV" ) ) != -1 )
        switch( c )
           {
            case  'a':   if ( strlen( optarg ) > MAX_ANCH_LEN )
                            {
                             fprintf( stderr, 
                                      "anchor enzyme length exceeds max of %d\n",
                                      MAX_ANCH_LEN );
                             fprintf( stderr,
                            "(consider recompiling with a larger valueof MAX_ANCH_LEN)\n");
                             exit( 1 );
                            }
                         anch_enz = uc( strdup( optarg ) ); 
                         
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
            case  'h':   help();
                         exit( 1 );
                         break;
            case  'l':   linker = uc( strdup( optarg ) );
                         break;
            case  'n':   seq_num = atoi( optarg ) - 1;
                         break;
            case  'F':   fragnumbers = 1;
                         break;
            case  'S':   seqnames = 1;
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

#define MAX_CBUFF  ( (2 * MAX_TAG_EXT) + MAX_ANCH_LEN )
char   cbuff[MAX_CBUFF+1];

int    enz_pos;
int    cbuff_len;

void  init_cbuff( void )
   {
    cbuff_len = 2 * tag_ext + anch_enz_len;  /* could check against MAX_CBUFF*/
    if ( verbose ) printf( "initializing cbuff, length %d..\n", cbuff_len );
    memset( cbuff, 'X', cbuff_len );
    cbuff[cbuff_len] = '\0';
    pos = -( anch_enz_len + tag_ext );
    if ( verbose ) printf( "cbuff initialized\n" );
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

void  print_tags( char *ftag, char *rtag, int tag_ext )
   {
    static char ft[MAX_TAG_EXT];
    static char rt[MAX_TAG_EXT];

    strncpy( ft, ftag, tag_ext );
    strncpy( rt, rtag, tag_ext );
/*    printf( "%s\t%2d\tr\t%d\n", rt, seq_num, pos );
    printf( "%s\t%2d\tf\t%d\n", ft, seq_num, pos ); */
    if ( seqnames )
       {
        printf( "%s        %s        r        %d", rt, header, pos );
        if ( fragnumbers )  printf( "       %10d", frag_num++ );
        printf( "\n" );
        printf( "%s        %s        f        %d", ft, header, pos );
        if ( fragnumbers )  printf( "       %10d", frag_num );
        printf( "\n" );
       }
    else
       {
        printf( "%s        %2d        r        %d", rt, seq_num, pos );
        if ( fragnumbers )  printf( "       %10d", frag_num++ );
        printf( "\n" );
        printf( "%s        %2d        f        %d", ft, seq_num, pos );
        if ( fragnumbers )  printf( "       %10d", frag_num );
        printf( "\n" );
       }
   }


void  handle_anchor( void )
   {
    static char rtag[MAX_TAG_EXT];

    if ( verbose ) printf( "     handling anchor\n" );

    rc( cbuff, tag_ext, rtag );
    print_tags( cbuff + tag_ext + anch_enz_len, rtag, tag_ext );
   }


/* handle( c ) - handle one input character  */

void  handle( char c )
   {
    if ( isalpha( c ) )
       {                              /* if its a nucleotide,    */
        insert( c );                  /* then put it into cbuff  */
        if ( is_anchor_site() )       /* and check to see if its */
            handle_anchor();          /* our anchor enzyme seqs  */
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
    char *h;

    if ( verbose ) printf( "new sequence\n" );

    flush_last_seq(); /* flush out last bit of buffer */
    init_cbuff();     /* new sequence => starting with fresh cbuff and fresh */
                      /* state variables */
    frag_num++;
    h = header;
    while ( (c = fgetc( f )) != EOF && c != '\n' )
       {
        if ( verbose ) putchar( c );
        if ( seqnames ) *h++ = c;
       }
    if ( verbose ) printf( "  char is now [%d]\n", c );
    if ( seqnames )
       {
        *h = '\0';
        h = header;
        while ( *h && *h != ' ' && *h != '|' )
            h++;
        *h = '\0';
       }    
    seq_num++;
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

    init_cbuff();

    if ( verbose )
       {
        printf( "anchor:     %s\n", anch_enz );
        printf( "anchor len: %d\n", anch_enz_len );
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
   }
