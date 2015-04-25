/* Program:      seqdiff.c                                                */
/* Programmer:   Sean McCorkle                                            */
/* Language:     C                                                        */
/* Description:  Compare two dna sequence strings, showing differences    */
/* $Id: seqdiff.c,v 3.8 2007/11/13 19:31:01 mccorkle Exp mccorkle $       */
/**************************************************************************/

#include <errno.h>
#include <libgen.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <string.h>
#include "seqlib.h"



                         /************************/
                         /* Macros and Constants */
                         /************************/

#undef   MAX_SEQ_LENGTH
#define  MAX_SEQ_LENGTH  8000  /* must override size in seqlib.h */
#define  MAX_FILE_NAME    100

                                /***********/
                                /* Globals */
                                /***********/


int  show_alignment     = 0;     /* set true by -a */
int  show_diffs         = 0;     /* set true by -l */
int  handle_translation = 0;     /* set true by -t */

static char seqdiff_rcs_id[] = 
      "$Id: seqdiff.c,v 3.8 2007/11/13 19:31:01 mccorkle Exp mccorkle $";

                         /*************************/
                         /* Function Declarations */
                         /*************************/

void  getargs( int argc, char **argv, char *file1, char *file2, FILE **a_file);
void  usage( void );
void  getseq( char *file, int *n, char *seq, char *hdr );


main( int argc, char **argv )

   {
    static char   file1[MAX_FILE_NAME];
    static char   file2[MAX_FILE_NAME];
    FILE         *al_file;
    int           n, m;                  /* lengths of seq1 and seq2 */
    static char   seq1[MAX_SEQ_LENGTH];
    static char   seq2[MAX_SEQ_LENGTH];
    static char   hdr1[MAX_HDR_LENGTH];
    static char   hdr2[MAX_HDR_LENGTH];

    static char   top[MAX_SEQ_LENGTH], 
                  mid[MAX_SEQ_LENGTH], 
                  bot[MAX_SEQ_LENGTH];
    static char   trans_top[MAX_SEQ_LENGTH], 
                  trans_bot[MAX_SEQ_LENGTH];
    int           l;                     /* length of alignment arrays */
    ALIGN_STATS   al_stats;

#ifdef METROWERKS
    show_alignment = 0;
    show_diffs = 0;
    strcpy( file1, "BB-045-C02D-R-472a-11-Jd" );
    strcpy( file2, "045-C02-RV-472a.11" );
#else
    getargs( argc, argv, file1, file2, &al_file );
#endif

    getseq( file1, &m, seq1, hdr1 );
#ifdef DEBUG
    printf( "sequence 1:\n" );
    print_seq( stdout, seq1, 0, 0, 0, 0, 0, 0 );
#endif

    getseq( file2, &n, seq2, hdr2 );
#ifdef DEBUG
    printf( "sequence 2:\n" );
    print_seq( stdout, seq2, 0, 0, 0, 0, 0, 0 );
    printf( "m, n = %d, %d\n", m, n );
#endif

    lpa_align( seq1, m, seq2, n, MAX_SEQ_LENGTH, handle_translation,
               top, mid, bot, 
               trans_top, trans_bot, &l );

#ifdef DEBUG
    printf( "seqdiff: align l is %d\n", l );
#endif
    get_align_stats( top, mid, bot, l, m, n, &al_stats, 
                     handle_translation, trans_top, trans_bot );
    print_align_stats( &al_stats, handle_translation );
    printf( "seq1 len:                %4d\n", m );
    printf( "seq2 len:                %4d\n", n );

    if ( show_diffs )
        print_diffs( &al_stats, handle_translation, top, mid, bot, 
                     trans_top, trans_bot, l );
    if ( show_alignment )
       {
        /*print_alignment( al_file, basename(file1), basename(file2), 
                         top, mid, bot, l );*/
        fprintf( al_file, "\n" );
        print_alignment( al_file, hdr1, hdr2, handle_translation,
                         top, mid, bot, trans_top, trans_bot, l );
       }

    /* Loose ends: not closing al_file */
   }


void getargs( int argc, char **argv, char *file1, char *file2, FILE **al_file )

   {
    char         *c;

    progname = strdup( basename( *argv ) );
    while ( --argc > 0 && **(++argv) == '-' )
       {
        c = (*argv)+1;
        if ( ! *c )
            usage();
        while ( *c )
            switch( *c++ )
	       {
	        case 'a':  show_alignment = 1;  
                           *al_file = fopen_safely( *++argv, "w" );
                           argc--;
                           break;
	        case 'l':  show_diffs = 1;      
                           break;
	        case 't':  handle_translation = 1;      
                           break;
                default:  usage();
	       }
       }
    if ( argc != 2 )
       usage();
    strncpy( file1, *argv++, MAX_FILE_NAME );
    strncpy( file2, *argv++, MAX_FILE_NAME );
#ifdef DEBUG
    printf( "file1: %s\n", file1 );
    printf( "file2: %s\n", file2 );
    printf( "show_alignment: %d, show_diffs %d\n",show_alignment, show_diffs );
#endif
   }


void usage( void )

   {
    printf( "usage: %s [-alt] <seq-file 1> <seq-file 2>\n", progname );
    exit( 1 );
   }


void getseq( char *file, int *n, char *seq, char *hdr )

   {
    int          seqf;
    static char  raw_seq[MAX_SEQ_LENGTH];

    if ( ( seqf = open_seqf( file, STD_ALLOWED, STD_IGNORED ) < 0 ) )
       {
        fprintf( stderr, "%s: fatal; ", progname );
        perror( file );
        exit( errno );
       }
    
    if ( read_seq( seqf, MAX_HDR_LENGTH, hdr, MAX_SEQ_LENGTH, raw_seq ) < 0 )
       {
        fprintf( stderr, "%s: fatal; ", progname );
        exit( 1 );
       }
    uppercase( raw_seq, seq );
    *n = strlen( seq );
    close_seqf( seqf );
   }



