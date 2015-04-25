/* Program:      n-mers.c                                                    */
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
/*                                                                           */
/* $Id$ */
/*                                                                           */
/*****************************************************************************/

#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#ifndef SOLARIS        /* Sun's are really a pain in the neck sometimes */
#include <unistd.h>
#endif
#include "seqlib.h"


#define  MAX_FILE_LEN    256    /* max length of file name */
#define  MAX_SEQ_LEN   10000    /* 200000*/    /* max length of one sequence */
#define  AVG_SEQ_LEN    2000    /* 8000*/   /* average sequence length */
#define  MAX_SEQUENCES  2000    /* 1000*/   /* max number of sequences */
#define  MAX_HDR_LEN     256    /* max length of header */

#define  HDR_BUFF_LEN  ( MAX_SEQUENCES * ( MAX_HDR_LEN + 1 ) )
#define  SEQ_BUFF_LEN  ( MAX_SEQUENCES * ( AVG_SEQ_LEN + 1 ) )

#define  PROTEINS     "arndcqeghilkmfpstwayvARNDCQEGHILKMFPSTWAYV"
#define  IGNORED      " \t\n*"
#define  HUGE_SIZE        10000000

static char rcs_id[] = 
      "$Id$";

                                /***********/
                                /* Globals */
                                /***********/


int      verbose      = 0;    /* verbose report for >thresh matches, -v opt */
int      length       = 18;   /* set by -l */

char     filename[MAX_FILE_LEN];

char     hdr_buff[ HDR_BUFF_LEN ];
char    *hdr_buff_end = hdr_buff + HDR_BUFF_LEN;
char    *next_hdr = hdr_buff;

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


/* this is really ugly.  Fix this! */

char *store_safely( char *data, char **start, char *end )
   {
    int  n;
    char *loc;

    n = strlen( data );
    if ( ( n + 1 ) >= end - *start )
        bailout( "store_safely limit failure" );
    loc = *start;
    (*start) += n + 1;
    return( strncpy( loc, data, n ) );
   }


void  parse_args( int argc, char **argv )
   {
    int          c;
    extern char *optarg;
    extern int   optind;

    while ( (c = getopt(argc, argv, "cdDmpt:s:T:v") ) != -1 )
        switch (c)
           {
            case 's':
                       length = atoi( optarg );
                       break;
            case 'v':
                       verbose = 1;
                       break;
           }

#ifdef DEBUG
    dump_sarray = 1;                 /* always on in DEBUG mode! */
    verbose = 1;                  
#endif
    strncpy( filename, argv[optind++], MAX_FILE_LEN );
   }



int  get_seq( int seqf, char *seq, int *len, char *hdr )
   {
    int  ret;

    ret = read_seq( seqf, MAX_HDR_LENGTH, hdr, MAX_SEQ_LEN, seq );
    if ( ret > 0 ) 
       {
        *len = strlen( seq );
        uppercase( seq, seq );
       }
    return( ret );
   }


void read_sequences ( char *filename )

   {
    int          seqf;
    static char  hdr[MAX_HDR_LENGTH];
    static char  seq[MAX_SEQ_LEN+1];
    int          seq_len;
    int          total_len = 0;
    int          i;

    char *allowed;
    char *ignored;

    if ( dna )
       {
        allowed = STD_ALLOWED;
        ignored = STD_IGNORED;
       }
    else      
       {
        allowed = PROTEINS;
        ignored = IGNORED;
       }

    if ( ( seqf = open_seqf( filename, allowed, ignored ) < 0 ) )
       {
        fprintf( stderr, "%s: fatal; ", progname );
        perror( filename );
        exit( errno );
       }
    while ( get_seq( seqf, seq, &seq_len, hdr ) > 0 )
       {
        if ( num_seqs >= MAX_SEQUENCES )
             bailout( "Exceeded maximum number of sequences" );
        header[num_seqs] = store_safely( hdr, &next_hdr, hdr_buff_end );
        sequence[num_seqs] = store_safely( seq, &next_seq, seq_buff_end );
        total_len += seq_len;
        num_seqs++;  
       }
    if ( verbose )
        printf( "%d sequences, total length is %d\n", num_seqs, total_len );
#ifdef DEBUG
    for ( i = 0; i < num_seqs; i++ )
       {
        printf( ">%s: length %d\n", header[i], strlen( sequence[i] ) );
        print_seq( stdout, sequence[i], 0, 0, 0, 0, 0, 0 );
       }
#endif
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

/* replace with binary search */

void lookup( int pos, int *id, int *offset )
   {
    int j = 0;

    while ( (j < num_seqs) & ((int) ( sequence[j] - seq_buff ) <= pos  ) )
       { 
        j++;
       }
    *id = --j;
    *offset = pos - (int) ( sequence[j] - seq_buff );
   }

int compare( char *a, char *b )
   {
    int l = 0;
    while ( *a && ( *a++ == *b++ ) )
        l++;
    return( l );
   }

#define MAX_PBUFF 60

#ifdef NADA_NAGA
void report( void )
   {
    int i, j;
    int l;
    int id1, id2;
    int offset1, offset2;
    int previous_above = 0;        /* 1 means previous above threshold too */
    int max_l = HUGE_SIZE;         /* max exact match throughout cluster */
                                   /* to show clusters */
    double      sum_pos = 0.0;
    double      sum_pos_sq = 0.0;
    double      avg_pos, sigma_pos;
    static char pbuff[MAX_PBUFF];
    static char fmt[256];

    if ( verbose )
        sprintf( fmt, "%%-%d.%ds [%%.%ds]\n", 
                      hdr_trunc, hdr_trunc, suffix_trunc );


    for ( i = 0, j = 1;  j < npos;   i++, j++ )
       {
        l = compare( seq_buff + pos[i], seq_buff + pos[j] );
        if ( l >= thresh )
           {
            if ( l < max_l )
                max_l = l;
    
/*            printf( "hit i,j = %d,%d l = %d\n%s %s\n", i, j, l,
	      seq_buff + pos[i], seq_buff + pos[j] ); */
/*            printf( "[%5d %5d   %5d %5d] ", i, pos[i], j, pos[j] );*/
            lookup( pos[i], &id1, &offset1 );
            lookup( pos[j], &id2, &offset2 );

            if ( print_matches )
               {
                printf( "%5d %5d (%5d) -- %5d (%5d)", 
                         l, id1, offset1, id2, offset2 );
                if ( previous_above > 0 )
                    printf( "      %4d", previous_above );
                else
                    printf( "          " );
                strncpy( pbuff, seq_buff + pos[j], 20 );
                pbuff[20] = '\0';
                printf( "  %20s", pbuff );
                printf( "\n" );
               }
            if ( clust_stats )
               {
                sum_pos += (double) offset1;
                sum_pos_sq += ((double) offset1) * ((double) offset1);
               }

            if ( verbose )
               {
                if ( ! previous_above )
                    printf( fmt, header[id1], seq_buff + pos[i] );
                printf( fmt, header[id2], seq_buff + pos[j] );
               }

            previous_above++;
           }
        else
           {
            if ( print_clust && ( previous_above >= clust_size ) )
               {
                strncpy( pbuff, seq_buff + pos[i], 
                            ( max_l > MAX_PBUFF ) ? MAX_PBUFF : max_l );
                pbuff[( max_l > MAX_PBUFF ) ? MAX_PBUFF : max_l ] = '\0';
                printf( "%8d %8d ", previous_above, max_l ); 
                if ( clust_stats )
                    {
                     /* tricky: using offset2 here*/
                     sum_pos += (double) offset2;
                     sum_pos_sq += (double) (offset2 * offset2);

                     avg_pos = sum_pos / (previous_above + 1);
                     sigma_pos = sqrt( (sum_pos_sq - (previous_above+1) * avg_pos 
                                        * avg_pos) / (previous_above) );

                     printf( "%6.1f %6.1f ", avg_pos, sigma_pos ); 
                    }
                printf( "%s\n", pbuff );
               }
            if ( print_clust && ( previous_above >= 0 ) )
                sum_pos = sum_pos_sq = 0.0;

            previous_above = 0;
            max_l = HUGE_SIZE;
           }
       }
   }
#endif


main( int argc, char **argv )
   {
    parse_args( argc, argv );

    if ( verbose )
        printf( "file is [%s]\n", filename );

    read_sequences( filename );

    create_suffix_array();

    print_suffix_array();
    
    /* report(); */

   }

