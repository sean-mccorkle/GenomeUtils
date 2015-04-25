/* Program:      k-mer-directory                                             */
/* Programmer:   Sean R. McCorkle                                            */
/* Language:     C                                                           */
/*                                                                           */
/* Description:  makes an index of k-mers and their positions within         */
/*               fasta sequences                                             */
/*                                                                           */
/* Usage:        k-mer-directory [-k <k>] <fasta-sequences>                  */
/*                                                                           */
/* Options:         -k <k>  k-mer size, default 21                           */
/*                                                                           */
/* $Id: k-mer-directory.c,v 0.1 2006/06/16 17:54:53 mccorkle Exp mccorkle $ */
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


#define  MAX_FILE_LEN     256   /* max length of file name */
#define  MAX_SEQ_LEN    20000   /* 200000*/    /* max length of one sequence */
#define  AVG_SEQ_LEN     6000   /* 8000*/   /* average sequence length */
#define  MAX_SEQUENCES    800   /* 1000*/   /* max number of sequences */
#define  MAX_HDR_LEN      256   /* max length of header */
#define  MAX_K            100   /* max value of k */

#define  HDR_BUFF_LEN  ( MAX_SEQUENCES * ( MAX_HDR_LEN + 1 ) )
#define  SEQ_BUFF_LEN  ( MAX_SEQUENCES * ( AVG_SEQ_LEN + 1 ) )


static char rcs_id[] = 
      "$Id: k-mer-directory.c,v 0.1 2006/06/16 17:54:53 mccorkle Exp mccorkle $";

                                /***********/
                                /* Globals */
                                /***********/

         /* options controlling operation */

int      k       = 21;        /* k-mer size, changed by -k<k> */

         /* options controlling output (debug or otherwise) */


int      verbose      = 0;    /* verbose report for >thresh matches, -v opt */
int      dump_sarray  = 0;    /* dump suffix array, set by -D option        */
int      suffix_trunc = 60;   /* truncate suffixes to this len. -T<n> option*/
int      hdr_trunc    = 20;

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
int      lengths[ MAX_SEQUENCES ];    /* one length for each sequence */
int      max_matched[ MAX_SEQUENCES ];  /* max matched so far */
char    *hdr_matched[ MAX_SEQUENCES ];
int      db_pos_matched[ MAX_SEQUENCES ];
int      sq_pos_matched[ MAX_SEQUENCES ];
int      num_seqs = 0;
int      max_db_seq_len = 0;

                               /****************/ 
                               /* suffix array */
                               /****************/

int      pos[ SEQ_BUFF_LEN ];    /* index into the pos array             */
int      npos;                   /* index of seq_buff                    */


void  bailout( char *msg )
   {
    fprintf( stderr, "%s: %s\n", progname, msg );
    exit( 1 );
   }

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

    while ( (c = getopt(argc, argv, "k:v") ) != -1 )
        switch (c)
           {
            case 'k':  k = atoi( optarg );
                       if ( k < 1 || k > MAX_K )
                          {
                           fprintf( stderr, "k must be between 1 and %d\n",
                                             MAX_K );
                           exit( 1 );
                          }
                       break;
            case 'v':
                       verbose = 1;
                       break;
           }

    strncpy( filename, argv[optind++], MAX_FILE_LEN );
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

int  get_seq( int seqf, char *seq, int *len, char *hdr )
   {
    int  ret;

    ret = read_seq( seqf, MAX_HDR_LEN, hdr, MAX_SEQ_LEN, seq );
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
    static char  hdr[MAX_HDR_LEN+1];
    static char  fhdr[MAX_HDR_LEN+1];
    static char  rhdr[MAX_HDR_LEN+1];
    static char  seq[MAX_SEQ_LEN+1];
    static char  rseq[MAX_SEQ_LEN+1];
    int          seq_len;
    int          total_len = 0;
    int          i;

    if ( verbose )
        printf( "read_sequences %s\n", filename );
    if ( ( seqf = open_seqf( filename, STD_ALLOWED, STD_IGNORED ) < 0 ) )
       {
        fprintf( stderr, "%s: fatal; ", progname );
        perror( filename );
        exit( errno );
       }
    while ( get_seq( seqf, seq, &seq_len, hdr ) > 0 )
       {
        if ( num_seqs >= MAX_SEQUENCES )
             bailout( "Exceeded maximum number of sequences" );
        fhdr[0] = 'f';
        fhdr[1] = ' ';
        strcpy( fhdr+2, hdr );
        header[num_seqs] = store_safely( fhdr, &next_hdr, hdr_buff_end );
        sequence[num_seqs] = store_safely( seq, &next_seq, seq_buff_end );
        total_len += seq_len;
        lengths[num_seqs] = seq_len;
        num_seqs++;

        /* reverse complement */
        if ( num_seqs >= MAX_SEQUENCES )
             bailout( "Exceeded maximum number of sequences" );
        dna_reverse_comp( seq, rseq );
        rhdr[0] = 'r';
        rhdr[1] = ' ';
        strcpy( rhdr+2, hdr );
        header[num_seqs] = store_safely( rhdr, &next_hdr, hdr_buff_end );
        sequence[num_seqs] = store_safely( rseq, &next_seq, seq_buff_end );
        total_len += seq_len;
        lengths[num_seqs] = seq_len;
        num_seqs++;

        if ( seq_len > max_db_seq_len )
            max_db_seq_len = seq_len; 
       }
    close_seqf( seqf );
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

void  init_best_matches ( void )
   {
    int id;

    for ( id = 0; id < num_seqs; id++ )
       {
        max_matched[id] = 0;
        db_pos_matched[id] = 0;
        sq_pos_matched[id] = 0;
        if ( !( hdr_matched[id] = (char *) malloc( MAX_HDR_LEN+1 ) ) )
           {
            fprintf( stderr, "out of memory allocating hdr_matched[%d]\n", id);
            exit(1);
           }
       }
   }

void  print_best_matches( void )
   {
    int id;

    for ( id = 0; id < num_seqs; id++ )
        if ( max_matched[id] > 0 )
           {
            printf( "%-20.20s %-20.20s %5d %10d %10d\n", 
                    header[id], hdr_matched[id], max_matched[id], 
                    db_pos_matched[id], sq_pos_matched[id] );
           }
   }

/* replace with binary search */

void id_lookup( int pos, int *id, int *offset )
   {
    int j = 0;

    while ( (j < num_seqs) & ((int) ( sequence[j] - seq_buff ) <= pos  ) )
       { 
        j++;
       }
    *id = --j;
    *offset = pos - (int) ( sequence[j] - seq_buff );
   }

void print_suffix_array( void )
   {
    int         i;
    static char fmt[256];
    int         id;
    int         offset;

    sprintf( fmt, "%%6d:  %%9d  [%%.%ds] %%8d %%.%ds\n", 
              suffix_trunc, hdr_trunc );

    for ( i = 0; i < npos; i++ )
       {
        id_lookup( pos[i], &id, &offset );
        printf( fmt, i, pos[i], seq_buff + pos[i], offset, header[id] );
       }
    printf( "done\n" );
   }


void  print_directory()
   {
    static char curr_kmer[MAX_K+1];
    int         curr_len = 0;
    int         i;
    int         id;
    int         offset;

    memset( curr_kmer, k+1, 0 );
    for ( i = 0; i < npos; i++ )
       {
        if ( strncmp( curr_kmer, seq_buff + pos[i], k ) != 0 )
           {
            strncpy( curr_kmer, seq_buff + pos[i], k );
            curr_len = strlen( curr_kmer );
            if ( curr_len == k )
                printf( "\n%s ", curr_kmer );
           }
        id_lookup( pos[i], &id, &offset );
        if ( curr_len == k )
            printf( "%4d ", offset );
       }
    printf( "\n" );
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


int compare( char *a, char *b )
   {
    int l = 0;
    while ( *a && ( *a++ == *b++ ) )
        l++;
    return( l );
   }



int sa_lookup( char *s, int thresh, int *m, int *id, int *offset )
   {
    int  u = 0;
    int  v = npos - 1;
    int  k;              /*   u < k < v */
    int  d;
    int  l, n;

    k =  ( u + v ) / 2;
    d = strcmp( s, seq_buff + pos[k] ); 
    while ( d != 0 && u <= v )              /* adjust k, u,v to match s[i] */
       {
        if ( d > 0 )
            u = k+1;  /* s > pos[k] */
        else
            v = k-1;   /* s < pos[k] */
        k =  ( u + v ) / 2;
        d = strcmp( s, seq_buff + pos[k] ); 
       }
    l  = compare( s, seq_buff + pos[k-1] );
    *m = compare( s, seq_buff + pos[k] );
    n  = compare( s, seq_buff + pos[k+1] );

    if ( l > *m )
       {
        *m = l;
        k = k-1;
       }
    if ( n > *m )
       {
        *m = n;
        k = k+1;
       }
    if ( *m >= thresh )
       {
        id_lookup( pos[k], id, offset );
        if ( header[*id][0] == 'r' )
           *offset = lengths[*id] - *offset;
        return( 1 );
       }
    else
        return( 0 );
   }

#ifdef NO

                       /**********************/
                       /* comparision buffer */
                       /**********************/

static char  cbuff[MAX_SEQ_LEN+1];
int          cbuff_len;

void  init_buff( void )
   {
    size_t i;

    cbuff_len = max_db_seq_len;
    if ( cbuff_len > MAX_SEQ_LEN )
       {
        fprintf( stderr, "sa_search: compare buffer length %d exceeds max %d\n",
                          cbuff_len, MAX_SEQ_LEN );
        exit( 1 );
       }
    
    for ( i = 0; i < cbuff_len; i++ )
        cbuff[i] = 'X';
    cbuff[cbuff_len] = '\0';

    if ( verbose )
        printf( "buffer initialized: %d\n", cbuff_len );
   }


void  enter_base( int c )
   {
    size_t i;

    memmove( cbuff, cbuff+1, cbuff_len-1 );  /* left shift forward tag buffer */

    cbuff[cbuff_len-1] = toupper( c );
   }

/*****************************************************/
/* scan_file() - open file and process its sequences */
/*****************************************************/

void  scan_file( char *filename )
   {
    FILE        *f;
    int          c, last;
    static char  curr_hdr[MAX_HDR_LEN+1];
    int          pos;   /* position within string */
    int          match_len;
    int          id, max_id;
    int          offset;

    if ( verbose )
        printf( "file %s\n", filename );
    f = open_file( filename );
    last = '\n';
    if ( verbose ) 
        printf( "max_db_seq_len is %d\n", max_db_seq_len );
    while ( (c = fgetc( f ) ) != EOF )
       {
        if ( last == '\n' && c == '>' )
           {
            fgets( curr_hdr, MAX_HDR_LEN, f );
            if ( verbose )
                printf( "seq: %s", curr_hdr );
            if ( strlen( curr_hdr ) > 0 )
                curr_hdr[strlen(curr_hdr)-1] = '\0';
            init_buff();
            /* pos = 1 - max_db_seq_len;  */ /* counting from one */
            pos = -max_db_seq_len;  /* counting from zero */
            last = '\n';
           }
        else if ( isalpha( c ) )
           {
            enter_base( c );
            ++pos;
            if ( verbose && (pos % 1000000 ) == 0 )
                printf( "%10d\n", pos );
            if ( sa_lookup( cbuff, thresh, &match_len, &id, &offset ) )
    	        if ( match_len > max_matched[id] )
                   {
                    if ( verbose )
                       printf( "id %d was %d now %d\n", 
                                    id, max_matched[id], match_len );
                    max_matched[id] = match_len;
                    db_pos_matched[id] = offset;
                    sq_pos_matched[id] = pos;
                    strncpy( hdr_matched[id], curr_hdr, MAX_HDR_LEN );
                   }
            last = c;
           }
        else
            last = c;
       }
    close_file( f );
   }
#endif

main( int argc, char **argv )
   {
    parse_args( argc, argv );

    if ( verbose )
        printf( "file is [%s]\n", filename );

    read_sequences( filename );

    create_suffix_array();

    print_directory();
   }

