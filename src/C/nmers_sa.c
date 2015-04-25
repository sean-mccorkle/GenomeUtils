/* Program:      nmers_sa.c                                                  */
/* Programmer:   Sean R. McCorkle                                            */
/* Language:     C                                                           */
/*                                                                           */
/* Description:  Lists n-mers and their frequencies in all input sequences   */
/*               (DNA or protein).                                           */
/*                                                                           */
/* Usage:        n-mers [-dtDTv] [ <file> [...]]                             */
/*                where file contains sequences (fasta format). "-" means    */
/*                scan std input.                                            */
/*                                                                           */
/* Options:                                                                  */
/*                                                                           */
/*                -l<n>  print n-mers of size <n>                            */
/*                -v     verbose output                                      */
/*                -D     dump suffix array (volumous output)                 */
/*                                                                           */
/* $Id$ */
/*                                                                           */
/*****************************************************************************/

#include <ctype.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#ifndef SOLARIS        /* Sun's are really a pain in the neck sometimes */
#include <unistd.h>
#endif


#define  MAX_FILE_LEN    256    /* max length of file name */
#define  AVG_SEQ_LEN  20000000    /* 8000*/   /* average sequence length */
#define  MAX_SEQUENCES    22    /* 1000*/   /* max number of sequences */
#define  MAX_HDR_LEN     256    /* max length of header */

/* #define  HDR_BUFF_LEN  ( MAX_SEQUENCES * ( MAX_HDR_LEN + 1 ) )*/
#define  SEQ_BUFF_LEN  ( MAX_SEQUENCES * ( AVG_SEQ_LEN + 1 ) )


int      verbose      = 0;    /* verbose report for >thresh matches, -v opt */
int      length       = 18;   /* set by -l */
int      dump_sarray  = 0;    /* set by -D */

char     filename[MAX_FILE_LEN];

/*
char     hdr_buff[ HDR_BUFF_LEN ];
char    *hdr_buff_end = hdr_buff + HDR_BUFF_LEN;
char    *next_hdr = hdr_buff;
*/

         /* every input string is stored in this array */

char     seq_buff[ SEQ_BUFF_LEN ];
char    *seq_buff_end = seq_buff + SEQ_BUFF_LEN;
char    *next_seq = seq_buff;

char    *header[ MAX_SEQUENCES ];     /* -> fasta header of sequence[i] */
char    *sequence[ MAX_SEQUENCES ];   /* -> each sequence[i]            */
int      num_seqs = 0;

                               /****************/ 
                               /* suffix array */
                               /****************/

int      pos[ SEQ_BUFF_LEN ];    /* index into the pos array             */
int      npos;                   /* index of seq_buff                    */

void  bailout( char *msg )
   {
    fprintf( stderr, "n-mers: %s\n", msg );
    exit( 1 );
   }

void  help( void )
   {
    fprintf( stderr, "What?\n" );
    exit( 1 );
   }

void  parse_args( int argc, char **argv, int *nfiles, char ***filenames )
   {
    int          c;
    extern char *optarg;
    extern int   optind;
    static char *def_files[] = { "-", "" };

    while ( (c = getopt(argc, argv, "Dl:v") ) != -1 )
        switch (c)
           {
            case 'D':  dump_sarray = 1;
                       break;
            case 'l':  length = atoi( optarg );
                       break;
            case 'v':  verbose = 1;
                       break;
            default:   help();
           }

#ifdef DEBUG
    dump_sarray = 1;                 /* always on in DEBUG mode! */
    verbose = 1;                  
#endif

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


void  load_sequences( char *filename )
   {
    FILE        *f;
    int          c;
    char         last_char;
    static char  hdrbuff[MAX_HDR_LEN+1];
    int          i;

    if ( ! (f = fopen( filename, "r" ) ) )
       {
        perror( filename );
        exit( errno );
       }
    if ( verbose )
        printf( "....loading %s\n", filename );

    last_char = '\n';
    while ( (c = fgetc( f )) != EOF )
       {
        if ( c == '>' && last_char == '\n' )   /* header line */
           {
            if ( num_seqs > 0 )             /* previous sequence? */
               {                            /* then close it off */
                *next_seq++ = '\0';
               }
            i = 0;
            while ( (c = fgetc( f )) != EOF  && c != '\n' && i < MAX_HDR_LEN)
                hdrbuff[i++] = c;
            if ( i >= MAX_HDR_LEN )
               {
                fclose( f );
                fprintf( stderr, 
                         "n-mers: (file %s) FASTA header exceeds %d\n",
                         MAX_HDR_LEN );
                bailout( "exiting" );
               }
            hdrbuff[i] = '\0';
            if ( ! (header[num_seqs] = strdup( hdrbuff )) )
                bailout( "unable to allocate header" );
            sequence[num_seqs] = next_seq;
            num_seqs++;
            if ( verbose )
                printf( "%s\n", hdrbuff );
            last_char = '\n';
           }
        else                                 /* sequence char */ 
           {
	    if ( isalpha( c ) )
               {
                if ( next_seq >= seq_buff_end )
                    bailout( "out of space for sequences" );
                *next_seq++ = c;
               }
            last_char = c;
           }
       }
    fclose( f );
   }



void print_suffix_array()
   {
    int         i;
    static char fmt[256];

    sprintf( fmt, "%%6d:  %%9d  [%%.%ds]\n", 60 );

    for ( i = 0; i < npos; i++ )
        printf( fmt, i, pos[i], seq_buff + pos[i] );
    printf( "done\n" );
   }



int  alpha_cmp( const void *a, const void *b )
   {
    return( strcmp( seq_buff + (*(int*)a), seq_buff +  (*(int*)b) ) );
   }

void  create_suffix_array()
   {
    int j;  /* position index of arrays */ 
    int i;  /* counts strings */

    i = j = npos = 0;
    if ( verbose )
        printf( "create suffix array\n" );
    while ( i < num_seqs )
       {
        if ( seq_buff[j] == '\0' )
           {
            i++;
#ifdef DEBUG
            printf( "create_suffix_array - string id now %d\n", i );
#endif
           }
        else 
            pos[npos++] = j;
        j++;
       }

#ifdef DEBUG
    print_suffix_array();
#endif

    if ( verbose )
        printf( "Starting qsort...\n" );
    qsort( pos, npos, sizeof( int ), alpha_cmp );
    if ( verbose )
        printf( "...Finished qsort\n" );

    if ( dump_sarray )
        print_suffix_array();
   }


void  report_n_mers( void )
   {
    int         i;
    int         k;
    int         n = npos - 1;
    static char fmt[256];

    sprintf( fmt, "%%4d %%.%ds\n", length );

    k = 1;
    i = 0;
    while ( i < n )
       {
        if ( strlen( seq_buff+pos[i] ) >= length )  /* improve this */
           {
            if ( strncmp( seq_buff+pos[i], seq_buff+pos[i+1], length ) == 0 )
                k++;
            else
               {
                printf( fmt, k, seq_buff+pos[i] );
                k = 1;
               }
           }
        else
           k = 1;
        i++;
       }
    printf( fmt, k, seq_buff+pos[i] );
   }


main( int argc, char **argv )
   {
    int     i;
    int     nfiles;
    char  **filenames;

    parse_args( argc, argv, &nfiles, &filenames );
    for ( i = 0; i < nfiles; i++ )
        load_sequences( filenames[i] );
    if ( verbose )
        for ( i = 0; i < num_seqs; i++ )
            printf( ">%s\n%s\n", header[i], sequence[i] );

    create_suffix_array();

    report_n_mers();
   }
