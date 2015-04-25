#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifndef SOLARIS        /* Sun's are really a pain in the neck sometimes */
#include <unistd.h>
#endif
#include "seqlib.h"

#define MAX_SEQ_LEN 5000   /* we're dealing mostly with single reads here */

/* globals which are potentially set by command line options */

double    thresh   = 2.0;       /* -p - 2% by default */
int       show_seq = 0;          /* -s */
int       n_files;
char    **files;

char      seq[MAX_SEQ_LEN];
int       seqlen;
char      hdr[MAX_HDR_LENGTH];
int       pos[MAX_SEQ_LEN];      /* positions of n's, plus end points */
int       npos;                  /* number of entries in pos */


/*******************************************************************************/
/* The usual option/argument parsing routine.  Set option globals accordingly, */
/* and put file arguments into global **files & set nfiles                     */
/*******************************************************************************/

void  parse_args( int argc, char **argv )
   {
    int          c;
    int          i;
    extern char *optarg;
    extern int   optind;
          
    while ( (c = getopt(argc, argv, "p:hs") ) != -1 )
        switch (c) 
           {
            case 'h': 
                       fprintf( stderr, "no help yet. sorry.\n" );
                       exit( 0 ); 
                       break;
            case 'p':
                       thresh = atof( optarg );
                       break;
            case 's': 
                       show_seq = 1;
                       break;
            case '?': 
                       fprintf( stderr, "Goodbye, Mr. Anderson\n" );
                       exit( 0 );
           }
    n_files = argc - optind;    
    files = (char **) malloc( n_files * sizeof( char * ) );
    i = 0;
    while ( optind < argc )
        files[i++] = strdup( argv[optind++] ); 
   }


/********************************************************/
/* This loads one sequence into the global seq[] buffer */
/********************************************************/

int  get_seq( int seqf )
   {
    int   ret;
    char *b;

    ret = read_seq( seqf, MAX_HDR_LENGTH, hdr, MAX_SEQ_LEN, seq );
    if ( ret > 0 ) 
       {
        seqlen = strlen( seq );
#ifdef DEBUG
        if ( b = index( hdr, ' ' ) )
  	    *b = '\0';                 /* cut off name at first blank */
#endif
       }
    return( ret );
   }


/**************************************************************************/
/* This sequential finds the positions of all the n's in the seq[] buffer */
/* and stores them in the array pos[].  In addition, the first element of */
/* pos[] is -1 (indicating the beginning of the sequence, and the last    */
/* element of pos[] is seqlen, indicating the end of the sequence         */
/**************************************************************************/

void  get_n_pos( void )
   {
    int  i;
    char *c;

    pos[0] = -1;
    npos = 1;
    for ( i = 0, c = seq;  i < seqlen;  i++, c++ )
        if ( *c == 'N' || *c == 'n' )
            pos[npos++] = i;

    pos[npos++] = seqlen;

#ifdef DEBUG
    printf( "Heres the pos array:\n" );
    for ( i = 0; i < npos; i++ )
        printf( "%d: %d\n", i, pos[i] );
#endif
   }


void  find_longest( 
                  int    *start,     /* left end of max length subsequence     */
                  int    *length,    /* length of max length subsequence       */
                  double *percnt     /* percentage N's of max len subsequence  */
                 )
   {
    int    i, j;          /* pos[i], pos[j] are left and right ends of  */
                          /*         candidate subsequence              */
    int    max_i;         /* left end of max candidate subsequence so far*/
    int    l;             /* length of current candidate */
    int    ncount;        /* number of N's in current candidate */
    double r;             /* percentage N's of current candidate */

    *start = -1;
    *length = 0;
    *percnt = -1.0;
    max_i = -1;

    for (  j = npos-1;  j > 0;  j--  )
        for (  i = 0;  i < j;  i++  )
           {
            l = pos[j] - ( pos[i] + 1 );
            if ( l > *length || (l == *length && i < *start) )
               {
                ncount = j - ( i + 1 );
                r = ( ncount * 100.0 ) / l; /* l must be > 0 if l > *length */
#ifdef DEBUG
                printf( "[%3d,%3d] %4d..%4d l = %4d, ncount %3d, %% = %5.2f\n",
                         i, j, pos[i], pos[j], l, ncount, r );
#endif
                if ( r <= thresh )
                   {
#ifdef DEBUG
                    printf( "*****biggest so far: l = %d\n", l );
#endif
                    *length = l;
                    max_i = i;
                    *percnt = r;
                   }
               }
           }
    if ( *length > 0 )
        *start = pos[max_i] + 1;
   }


/***************************************************************/
/* Output the subsequence of seq[] defined by start and length */
/***************************************************************/

void  output_cut( int start, int length, double percnt )
   {
    printf( ">%s  start: %d  length: %d  %5.2lf%% N's\n", 
              hdr,    start,     length,    percnt         );
    if ( length > 0 )
        print_seq( stdout, seq, start, start + length, 0, 0, 0, 0 );
   }


/*****************************************************************************/
/* This handles everything for the sequence current in the global seq buffer */
/*****************************************************************************/

void  process_seq( void )
   {
    int    length, start;
    double percnt;

#ifdef DEBUG
    printf( "Processing seq [%d] %s\n", seqlen, hdr );
#endif

    get_n_pos();
    find_longest( &start, &length, &percnt );
    output_cut( start, length, percnt );
   }

 
                                /****************/                            
                                /* Main Program */
                                /****************/                            

int  main( int argc, char **argv )
   {
    int  i;
    int  seqf;

    parse_args( argc, argv );

    for ( i = 0; i < n_files; i++ )
       {
        if ( ( seqf = open_seqf( files[i], STD_ALLOWED, STD_IGNORED ) < 0 ) )
           {
            fprintf( stderr, "%s: fatal; ", progname );
            perror( files[i] );
            exit( errno );
           }
        while ( get_seq( seqf ) > 0 )
            process_seq();
        close_seqf( seqf );
       }
   }



