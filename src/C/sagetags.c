/* Program:     sagetags.c                                                   */
/* Programmer:  Sean R McCorkle                                              */
/*              Biology Department, Brookhaven National Laboratory           */
/* Language:    C                                                            */
/*                                                                           */
/* Description: extracts SAGE tag sequences from mRNA transcripts            */
/*                                                                           */
/* Usage:       sagetags [-aehlrs] [file ...]                                */
/*                                                                           */
/*              where [files] are mRNA/cDNA transcripts, in FASTA format.    */
/*              stdin is scanned it no filemames are specified, and the      */
/*              filename "-" can also be used to indicate stdin.             */
/*                                                                           */
/* Options:       -a       show ALL tags from ALL enzyme sites in each       */
/*                         transcript (in both directions)                   */
/*                -e<seq>  anchor restriction enzyme sequence (default CATG) */
/*                -h       print help, then exit                             */
/*                -l<n>    extend tags <n> bases after enzyme (default 10)   */
/*                -r       show tags in both forward and reverse directions  */
/*                -s       print complete suffixes, from last enzyme site to */
/*                         3' end of transcript                              */
/*                                                                           */
/* Output:      one line per tag, four columns each:                         */
/*                                                                           */
/*              tag-sequence  dir  pos  transcript-id                        */
/*                (the last is the fasta header)                             */
/*                                                                           */
/* Example output from  sagetags -r                                          */
/*                                                                           */
/*    seq        dir    pos transcript id                                    */
/*        .      .       .        .                                          */
/*        .      .       .        .                                          */
/*    AAGTTGCCTA f      923 NP040439 similar to Medicago truncatula  MtN21   */
/*    TTAGAAACAA r      966 NP040439 similar to Medicago truncatula  MtN21   */
/*    TCCATCAGGA f      972 NP040494 ERG protein                             */
/*    AAATTAGTCA r      702 NP040494 ERG protein                             */
/*    AGTCAAAAAA f     2795 NP040620 respiratory burst oxidase protein E     */
/*    AAGTTCATCT r     2072 NP040620 respiratory burst oxidase protein E     */
/*    CAAAGAACGG f     2376 NP040621 respiratory burst oxidase protein C     */
/*        .      .       .        .                                          */
/*        .      .       .        .                                          */
/*                                                                           */
/* Compiling:                                                                */
/*             cc -O -o sagetags sagetags.c                                  */
/*                                                                           */
/*             should do the trick on most unix systems (MacOSX, linux etc)  */
/*                                                                           */
/* $Id: sagetags.c,v 0.4 2005/02/27 19:41:50 mccorkle Exp mccorkle $         */
/*****************************************************************************/

#include <ctype.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static char stags_rcs_id[] = 
    "$Id: sagetags.c,v 0.4 2005/02/27 19:41:50 mccorkle Exp mccorkle $";

#define  MAX_TAG_LEN 40

/* globals which are set or modified by command line options in parse_args() */

int    all_tags       = 0;       /* set by -a */
char  *re_seq         = "CATG";  /* set by -e */
int    tag_extent     = 10;      /* set by -l */
int    reverse_tags   = 0;       /* set by -r */
int    print_suffixes = 0;       /* set by -s */
int    n_files;
char **files;

/* other globals */

#define  POLY_A_LEN  50          /* len. of poly_A tail; must be >MAX_TAG_LEN*/
char    *poly_a;                 /* we force a polyA tail on end of each mRNA*/
int      re_len;                 /* length of punctuating enzyme sequence    */
int      tag_len;                /* total length of tag, included enzyme seq */

#define  MAX_SEQ_LEN  200000
char     sequence[MAX_SEQ_LEN+1];  /* buffer to hold one transcript */
int      seq_len;

#define  MAX_HDR_LEN  255          /* max fasta header length */
char     header[MAX_HDR_LEN+1];

void  usage( void )     /* print help, then exit */
   {
    printf( " \n\
Usage:       sagetags [-aehlrs] [file ...]                                \n\
                                                                          \n\
             where [files] are mRNA/cDNA transcripts, in FASTA format.    \n\
             stdin is scanned it no filemames are specified, and the      \n\
             filename \"-\" can also be used to indicate stdin.           \n\
                                                                          \n\
Options:       -a       show ALL tags from ALL enzyme sites in each       \n\
                        transcript (in both directions)                   \n\
               -e<seq>  anchor restriction enzyme sequence (default CATG) \n\
               -h       print help, then exit                             \n\
               -l<n>    extend tags <n> bases after enzyme (default 10)   \n\
               -r       show tags in both forward and reverse directions  \n\
               -s       print complete suffixes, from last enzyme site to \n\
                        3' end of transcript                              \n\
                                                                          \n\
Output:      one line per tag, four columns each:                         \n\
                                                                          \n\
             tag-sequence  dir  pos  transcript-id                        \n\
               (the last is the fasta header)                             \n\
                                                                          \n\
Example output from  sagetags -r                                          \n\
                                                                          \n\
   seq        dir    pos transcript id                                    \n\
       .      .       .        .                                          \n\
       .      .       .        .                                          \n\
   AAGTTGCCTA f      923 NP040439 similar to Medicago truncatula  MtN21   \n\
   TTAGAAACAA r      966 NP040439 similar to Medicago truncatula  MtN21   \n\
   TCCATCAGGA f      972 NP040494 ERG protein                             \n\
   AAATTAGTCA r      702 NP040494 ERG protein                             \n\
       .      .       .        .                                          \n\
       .      .       .        .                                          \n\
                                                                          \n\
" );
    exit( 0 );
   }


/* uc( str ) - convert null-terminated str into uppercase, return ptr to s */

char *uc( char *s )
   {
    char *t;
    
    for ( t = s; *t; t++ )
        *t = toupper( *t );
    return( s );
   }


void  parse_args( int argc, char **argv )
   {
    int          c;
    int          i;
    extern char *optarg;
    extern int   optind;
    static char *def_files[] = { "-", "" };
          
    while ( (c = getopt(argc, argv, "ae:l:hrs") ) != -1 )
        switch (c) 
           {
            case 'a':
                       all_tags = 1;
                       break;
            case 'e':
                       re_seq = uc( strdup( optarg ) );
                       break;
            case 'l':
                       tag_extent = atoi( optarg );
                       if ( tag_extent < 0 || tag_extent > MAX_TAG_LEN )
                          {
                           fprintf( stderr, 
                              "tag extent %d must be between 0 and %d\n", 
                              tag_extent, MAX_TAG_LEN );
                           exit( 1 );
                          }
                       break;
            case 'h': 
                       usage();
                       break;
            case 'r': 
                       reverse_tags = 1;
                       break;
            case 's': 
                       print_suffixes = 1;
                       break;
            default:
                       fprintf( stderr, 
                                "bad option; (sagetags -h gives help)\n", c );
                       exit( 1 );
                       break;
           }
    argc -= optind;
    argv += optind;
    if ( argc > 0 )
       {
        n_files = argc;
        files = argv;
       }
    else
       {
        n_files = 1;
        files = def_files;
       }
   }



/* allocates a string and fills it with n copies of character c */

char *replicate( char c, int n )
   {
    char *s, *t;
    t = s = malloc( n+1 );
    while ( n-- > 0 )
        *t++ = c;
    *t = '\0';
    return( s );
   }



/* print out a tag hit */

void  print_tag( char *s, int tag_extent, char dir, int pos, char *name  )
   {
    static char    tag[MAX_TAG_LEN+1];

    strncpy( tag, s, tag_extent );
    printf( "%s %c %8d %s\n", tag, dir, pos, name );
   }

/* look for occurances of re_seq in string seq - report tag if found   */
/* if all_tags is set, all tags are report, otherwise only the last is */

void   tag_search( char *re_seq, int re_len, int tag_extent, 
                   char *seq, char dir, char *name )
   {
    char          *p;      /* position of punc. restriction site in seq */
    char          *last_p;

    p = seq;
    last_p = NULL;
    while ( p = strstr( p, re_seq ) )
       {
        if ( all_tags )
            print_tag( p+re_len, tag_extent, dir, (int)(p - seq), name );
        last_p = p++;
       }
    if ( last_p && ! all_tags )
        print_tag( last_p+re_len, tag_extent, dir, (int)(last_p - seq), name );
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


void clear_seq( void )           /* reset sequence buffer for new seq */
   {
    seq_len = 0;
   }


void append_seq( char c )        /* put c on end of sequence buffer */
   {
    if ( seq_len < MAX_SEQ_LEN )
        sequence[seq_len++] = c;
    else
       {
        fprintf( stderr, "mRNA/cDNA transcript length exceeds %d\n", 
                          MAX_SEQ_LEN );
        fprintf( stderr,"(consider increasing MAX_SEQ_LEN and recompiling)\n");
        exit( 1 );
       }
   }

/* scan sequence s for tags, in one direction */

void process_one_direction( char *s, int len, char dir )
   {
    static char seq[MAX_SEQ_LEN+POLY_A_LEN+1];

    strncpy( seq, s, MAX_SEQ_LEN );
    strncat( seq, poly_a, POLY_A_LEN );   /* append polyA tail here */
    tag_search( re_seq, re_len, tag_extent, seq, dir, header );
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


void rc_in_situ( char *seq )   /* reverse complement seq, leaving result in  */
   {                           /* same location                              */
    char *p, *q;
    char c;

    p = seq;
    q = seq + strlen( seq ) - 1;

    while ( q > p )
       {
        c = *p;
        *p++ = comp_char[*q];
        *q-- = comp_char[c];
       }        
        
    if ( q == p )
       *p = comp_char[*p];
   }


void process_seq( void )
   {
    if ( seq_len > 0 )
       {
        sequence[seq_len] = '\0';
        process_one_direction( sequence, seq_len, 'f' );
        if ( all_tags || reverse_tags )
           {
            rc_in_situ( sequence );
            process_one_direction( sequence, seq_len, 'r' );
           }
       }
   }

void get_fasta_header( FILE *f )    /* reads one line worth of data, */
   {                                /* puts into header[]            */
    int c;
    int  i = 0;

    c = fgetc( f );
    while ( c != EOF && c != '\n' )
       {
        if ( i < MAX_HDR_LEN )     /* ignore header past MAX_HDR_LEN */
            header[i++] = c;
        c = fgetc( f );
       }
    header[i] = '\0';
   }

/* read one file of transcript sequences, and scan for tags & print */

void  process_file ( char *filename )
   {
    FILE   *f;
    int     c;
    char    last_c;

    f = open_file( filename );
    last_c = '\n';
    while ( ( c = fgetc( f ) ) != EOF )
       {
        if ( c == '>' && last_c == '\n' )
           {
            process_seq();         /* handle previous sequence, if any */
            clear_seq();
            get_fasta_header( f );
           }
        else if ( isalpha( c ) )
           {
            append_seq( c );
            last_c = c;
           }
        else
            last_c = c;
       }
    close_file( f );    
    process_seq();    /* don't forget last seq in buffer!*/
   }


                              /****************/
                              /* Main Program */
                              /****************/


int  main( int argc, char **argv )
   {
    int            i;

    parse_args( argc, argv );               /* get filenames, options from */
                                            /* command line                */
    re_len = strlen( re_seq );
    tag_len = re_len + tag_extent;
    poly_a = replicate( 'A', POLY_A_LEN );  /* "AAAAAAAAAAAAAAAAAAAA..." */

    for ( i = 0; i < n_files; i++ )         /* handle each file on command */
        process_file( files[i] );           /* line                        */
   }



