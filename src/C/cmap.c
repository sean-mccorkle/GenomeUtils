#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "seqlib.h"



#define  MAX_FILE_LEN    256    /* max length of file name */
#define  MAX_SEQ_LEN  200000

char     filename1[MAX_FILE_LEN];
char     filename2[MAX_FILE_LEN];

char     seq[MAX_SEQ_LEN+1];
int      seq_len;
int      pos[MAX_SEQ_LEN];

char     seq2[MAX_SEQ_LEN+1];
int      seq2_len;


int  get_seq( int seqf, char *seq, int *len )
   {
    int            ret;
    static char    hdr[MAX_HDR_LENGTH];

    ret = read_seq( seqf, MAX_HDR_LENGTH, hdr, MAX_SEQ_LEN, seq );
    if ( ret > 0 ) 
       {
        *len = strlen( seq );
        uppercase( seq, seq );
       }
    return( ret );
   }

void  parse_args( int argc, char **argv )
   {
    strncpy( filename1, argv[1], MAX_FILE_LEN );
    strncpy( filename2, argv[2], MAX_FILE_LEN );
   }

void  read_2nd_set( char *filename2 )
   {
    int            seqf;

    if ( ( seqf = open_seqf( filename2, "acgtnxACGTNX", STD_IGNORED ) < 0 ) )
       {
        fprintf( stderr, "%s: fatal; ", progname );
        perror( filename2 );
        exit( errno );
       }
    if ( get_seq( seqf, seq2, &seq2_len ) <= 0 )
       {
        fprintf( stderr, "%s: no sequence in file %s\n", progname, filename2 );
        exit( 1 );
       }

    close_seqf( seqf );
   }

int  alpha_cmp( const void *a, const void *b )
   {
    return( strcmp( seq + (*(int*)a), seq +  (*(int*)b) ) );
   }

void  create_suffix_array( char *s, int len, int *pos )
   {
    int p;

    for ( p = 0; p < len; p++ )
        pos[p] = p;
    qsort( pos, len, sizeof( int ), alpha_cmp );
   }

void print_suffix_array( int *pos, int seq_len )
   {
    int i;

    for ( i = 0; i < seq_len; i++ )
        printf( "%d:  %d  [%s]\n", i, pos[i], seq + pos[i] );
#ifdef NO
        printf( "%d:  %d\n", i, strlen( seq + pos[i] ) );
#endif
    printf( "done\n" );
   }

int min( int a, int b )
   { 
    return( ( a < b ) ? a : b );
   }


void  search_suffix_array( int *pos, int seq_len, char *pat )
   {
    int l, r, m;
    int cmp, pl;

    pl = strlen( pat );    
    printf( "Searching for [%s] %d\n", pat, pl );
    l = 0;
    r = seq_len;
    while ( l < r )
        {
         m = ( r + l ) / 2;
         cmp = strncmp( pat, seq + pos[m], min( pl, seq_len - pos[m] ) );
         printf( "checking %d %d %d [%s]: %d\n", r, l, m, seq + pos[m], cmp );
         if ( cmp == 0 )
             l = r = m;
         else if ( cmp < 0 )
             r = m;
         else       
             l = m;
        }
    printf( "l = %d, pos %d, [%s] : [%s]\n", l, pos[l], pat, seq+pos[l] );
   }


int thresh_search( int *pos, int seq_len, char *pat, int thresh,
                   int *where, int *len )
   {
    int l, r, m;
    int cmp, pl;
    int l_len, r_len, mlr;
    int j, k;

    printf( "Searching for [%s] %d, thresh\n", pat, thresh );
    l = 0;
    r = seq_len;
    l_len = r_len = mlr = 0;
    while ( r - l >= 2 )
        {
         m = ( r + l ) / 2;
         printf( "....l = %d  r = %d  ->  m = %d\n", l, r, m );
         /*compare chars at m, starting at mlr */
         j = mlr;
         k = pos[m] + j;
         while ( pat[j] != '\0' && k < seq_len && pat[j] == seq[k] )
            { j++; k++; }
         /*either 1) pat ran out - found
                  2) seq ran out ?
                  3) neither ran out - one is lex before the other
	  */
         printf( "          j = %d  k = %d,  pat[j] = %c seq[k] = %c\n",
                              j, k, pat[j], seq[k] );
         if ( pat[j] == '\0' || k >= seq_len )
             l = r = m;
         else if ( pat[j] < seq[k] )
            {
             printf( "              less\n" );
             r = m;
             r_len = j;
            }
         else       
            {
             printf( "              not less\n" );
             l = m;
             l_len = j;
            }
         mlr = min( l_len, r_len );
        }

    printf( "l = %d, pos %d, [%s] : [%s]\n", l, pos[l], pat, seq+pos[l] );
    printf( "mlr = %d, m = %d, k = %d j = %d\n", mlr, m, k, j );
    if ( j >= thresh )
       {
        *where = m;
        *len = j;
        return( 1 );
       }
    else
        return( 0 );
   }


main( int argc, char **argv )
   {
    int            seqf;
    int            where, len;
    int            thresh = 5;
    char  *p;

    parse_args( argc, argv );

    printf( "file1 is [%s]\n", filename1 );
    printf( "file2 is [%s]\n", filename2 );

    read_2nd_set( filename2 );
    if ( ( seqf = open_seqf( filename1, "acgtnxACGTNX", STD_IGNORED ) < 0 ) )
       {
        fprintf( stderr, "%s: fatal; ", progname );
        perror( filename1 );
        exit( errno );
       }
    if ( get_seq( seqf, seq, &seq_len ) <= 0 )
       {
        fprintf( stderr, "%s: no sequence in file %s\n", progname, filename1 );
        exit( 1 );
       }
    printf( "length is %d\n", seq_len );

    print_seq( stdout, seq, 0, 0, 0, 0, 0, 0 );


    create_suffix_array( seq, seq_len, pos );

    print_suffix_array( pos, seq_len );

    for ( p = seq2; *p != 0; p++ )
        if ( thresh_search( pos, seq_len, p, thresh, &where, &len ) )
            printf( "find @ %d len %d\n", where, len );
        else
            printf( "No find\n" );

   }

