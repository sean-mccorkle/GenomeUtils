#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifndef SOLARIS        /* Sun's are really a pain in the neck sometimes */
#include <unistd.h>
#endif
#include "seqlib.h"


#define  MAX_SEQ_LEN  2000000
#define  MAX_TAG_SIZE 20

/* globals which are potentially set by command line options */

int    all_tags       = 0;       /* -a */
char  *re_seq         = "CATG";  /* -e */
int    tag_extent     = 10;      /* -l */
int    verbose        = 0;       /* -v */
int    print_suffixes = 0;       /* -s */
int    n_files;
char **files;
char  *notag;   /* "XXXXX...." */

/* globals for calculations */


#define MAX_SEQUENCES 30000

char *seqnames[MAX_SEQUENCES];



void  parse_args( int argc, char **argv )
   {
    int          c;
    int          i;
    extern char *optarg;
    extern int   optind;
          
    while ( (c = getopt(argc, argv, "ae:l:hstTv") ) != -1 )
        switch (c) 
           {
            case 'a':
                       all_tags = 1;
                       break;
            case 'e':
                       re_seq = strdup( optarg );
                       uppercase( re_seq, re_seq );
                       break;
            case 'l':
                       tag_extent = atoi( optarg );
                       break;
            case 'h': 
                       fprintf( stderr, "no help yet. sorry.\n" );
                       exit( 0 ); 
                       break;
            case 's': 
                       print_suffixes = 1;
                       break;
            case 'v': 
                       verbose = 1;
                       break;
            case '?': 
                       fprintf( stderr, "Goodbye, Mr. Anderson\n" );
                       exit( 0 );
           }
    n_files = argc - optind;    
    files = (char **) malloc( n_files * sizeof( char * ) );
    i = 0;
    while( optind < argc )
        files[i++] = strdup( argv[optind++] ); 
   }


int  get_seq( int seqf, char *seq, int *len, char *hdr )
   {
    int   ret;
    char *b;

    ret = read_seq( seqf, MAX_HDR_LENGTH, hdr, MAX_SEQ_LEN, seq );
    if ( ret > 0 ) 
       {
        *len = strlen( seq );
        uppercase( seq, seq );
#ifdef DEBUG
        if ( b = index( hdr, ' ' ) )
  	    *b = '\0';                 /* cut off name at first blank */
#endif
       }
    return( ret );
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

void  find_tag( char *re_seq, int re_len, int tag_extent,
                 char *seq, char *tag )

   {
    char          *p;      /* position of punc. restriction site in seq */
    char          *last_p;

    strcpy( tag, notag );
    p = seq;
    last_p = NULL;
    while ( p = strstr( p, re_seq ) )
        last_p = p++;
    if ( last_p )
        strncpy( tag, last_p + re_len, tag_extent );
    /*printf( "tag is [%s]\n", tag ); */
    
   }

void  find_all_tags( char *re_seq, int re_len, int tag_extent,
                     char *seq, char dir, char *name )

   {
    char          *p;      /* position of punc. restriction site in seq */
    static char    tag[MAX_TAG_SIZE];
    int            i;
    
    p = seq;
    i = 0;
    while ( p = strstr( p, re_seq ) )
       {
        strncpy( tag, p + re_len, tag_extent );
        if ( dir == 'f' )
            printf( "%s  %s  {%c%d} %s\n", tag, notag, dir, ++i, name );
        else
            printf( "%s  %s  {%c%d} %s\n", notag, tag, dir, ++i, name );
        p++;
       }
   }


int  main( int argc, char **argv )
   {
    int            seqf;
    int            i;
    static char    seq[MAX_SEQ_LEN];
    static char    rseq[MAX_SEQ_LEN];
    int            len;    /* length of sequence being processed */
    int            re_len; /* length of punctuating enzyme sequence */
    static char    ftag[200];
    static char    rtag[200];
    int            tag_len;
    static char    seqnam[MAX_HDR_LENGTH];
    static char   *poly_a;    /* "AAAAAAAAAAAAAAAAAAAA"; */

    poly_a = replicate( 'A', 50 );
    parse_args( argc, argv );
    re_len = strlen( re_seq );
    tag_len = re_len + tag_extent;
    notag = replicate( 'X', tag_extent );
    for ( i = 0; i < n_files; i++ )
       {
        if ( ( seqf = open_seqf( files[i], STD_ALLOWED, STD_IGNORED ) < 0 ) )
           {
            fprintf( stderr, "%s: fatal; ", progname );
            perror( files[i] );
            exit( errno );
           }
        while ( get_seq( seqf, seq, &len, seqnam ) > 0 )
           {
            dna_reverse_comp( seq, rseq );
            strcat( seq, poly_a ); /* fix this! */
            strcat( rseq, poly_a ); /* fix this! */
#ifdef NO
            printf( "sequence loaded\n" );
            print_seq( stdout, seq, 0, 0, 5, 0, 0, 0 );
            printf( "Rc\n" );
            print_seq( stdout, rseq, 0, 0, 5, 0, 0, 0 );
#endif
            if ( all_tags )
               {
                find_all_tags( re_seq, re_len, tag_extent, seq, 'f', seqnam );
                find_all_tags( re_seq, re_len, tag_extent, rseq, 'r', seqnam );
               }
            else
               {
                find_tag( re_seq, re_len, tag_extent, seq, ftag );
                find_tag( re_seq, re_len, tag_extent, rseq, rtag );
                printf( "%s  %s  %s\n", ftag, rtag, seqnam );
               }
           }
        close_seqf( seqf );
       }
   }





