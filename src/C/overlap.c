/* Program:      overlap.c                                                   */
/* Programmer:   Sean R. McCorkle                                            */
/* Language:     C                                                           */
/*                                                                           */
/* Description:  Read two sequences (fasta format) and print position of     */
/*               maximum overlap between the two.                            */
/*                                                                           */
/* Usage:        overlap [-hs] [-r <repeatfile>] <file1> <file2>             */
/*                                                                           */
/* Options:        -h       print help then exit                             */
/*                 -s       print sequences                                  */
/*                 -r<file> exclude sequences found in this file (repeats)   */
/*                                                                           */
/* Note to self: change lengths from size_t to whatever (find out what       */
/*   whatever is output offsets are starting from zero.  change?             */
/*                                                                           */
/* $Id: overlap.c,v 0.5 2003/12/04 22:42:57 mccorkle Exp mccorkle $ */
/*                                                                           */
/*****************************************************************************/

#include <ctype.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

static char overlap_rcs_id[] =
     "$Id: overlap.c,v 0.5 2003/12/04 22:42:57 mccorkle Exp mccorkle $";

#define MAX_FILENAME 256
#define MAX_SEQ_LEN  100000
#define MAX_HDR_LEN  256


int   print_seq = 0;     /* set by -s */
 
void  usage( void )
   {
    fprintf( stderr, "usage:  overlap [-hs] [-r<file>] <file1> <file2>\n" );
   }

void  parse_args( int argc, char **argv, char *file1, char *file2,
                  char *repeatfile )
   {
    int          c;
    extern char *optarg;
    extern int   optind;

    *repeatfile = '\0';
    while ( (c = getopt( argc, argv, "hr:s" )) != -1 )
        switch ( c )
           {
            case 'r': strncpy( repeatfile, optarg, MAX_FILENAME );
                      break;
            case 's': print_seq = 1;
                      break;
            case 'h': usage();
                      exit( 0 );
            default:  usage();
                      exit( 1 );
           }

    argc -= optind; 
    argv += optind;   
    if ( argc != 2 )
       {
        usage();
        exit( 1 );
       }
    strncpy( file1, argv[0], MAX_FILENAME );
    strncpy( file2, argv[1], MAX_FILENAME );
   }


int  filesize( char *file )             /* returns the size of file in bytes*/
   {
    struct stat s;

    if ( stat( file, &s ) != 0 )
       {
        perror( file );
        exit( errno );
       }
    return( s.st_size );
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




/* c = get_seq( f, buff, max );   */
/* if max is exceeded, error exit */

int get_seq( FILE *f, char *buff, int max )  
   {
    char  *p;
    char  *p_end;
    int    c;

    p = buff;
    p_end = p + max;
    while ( (c = fgetc( f )) != EOF && c != '>' && p < p_end )
       if ( isalpha( c ) )
           *p++ = toupper( c );
    if ( p >= p_end )
       {
        fprintf( stderr, "get_seq: buffer overrun\n" );
        exit( 1 );
       }
    *p++ = '\0';
    return( c );
   }

  /* read fasta sequence from  */
 /* filename and put into buff*/

int  load_seq( char *filename, char *buff, int max ) 
   {                                           
    int          c;
    FILE        *f;
    static char  hdr[MAX_HDR_LEN+1];
   
    if ( !(f = fopen( filename, "r" )) )
       {
        perror( filename );
        exit( errno );
       }
    fgets( hdr, MAX_HDR_LEN, f );
    if ( hdr[0] != '>' )
       {
        fprintf( stderr, "%s: first line is not FASTA hdr\n", filename );
        fclose( f );
        exit( 1 );
       }
    get_seq( f, buff, max );
    fclose( f );
    return( strlen( buff ) );
   }


#define  MAX_REPEATS 100
#define  MAX_REP_LEN 2000

int    n_repeats = 0;
char  *repeats[MAX_REPEATS];

/* note to self: this should really use some common routine with load_seq */

void  load_repeats( char *file )
   {
    FILE        *f;
    int          c;
    static char  buff[MAX_REP_LEN];
    static char  hdr[MAX_HDR_LEN+1];  /* not used; maybe later */

    /* printf( "Loading repeats from [%s]\n", file ); */

    if ( !(f = fopen( file, "r")) )
       {
        perror( file );
        exit( errno );
       }
    c = fgetc( f );
    n_repeats = 0;
    while ( n_repeats < MAX_REPEATS && c == '>' )
       {
        fgets( hdr, MAX_HDR_LEN, f );        /* discard fasta hdr for now */
        c = get_seq( f, buff, MAX_REP_LEN );
        if ( !( repeats[n_repeats++] = strdup( buff ) ) )
           {
            perror( "failure to allocate repeat seq storage" );
            exit( errno );
           }
        if ( !( repeats[n_repeats++] = strdup( buff ) ) )
           {
            perror( "failure to allocate repeat seq storage" );
            exit( errno );
           }
        rc_in_situ( repeats[n_repeats-1] );
       }
    
    if  ( n_repeats >= MAX_REPEATS )
       {
        fprintf( stderr, "too many repeat sequences\n" );
        exit( 1 );
       }
    fclose( f );
   }


int  alpha_cmp( const void *a, const void *b )     /* string comparision for */
   {                                               /* use in qsort()         */
    return( strcmp( *(char **)a, *(char **)b ) );
   }


create_suffix_array( char *sa[], char *seqbuff, size_t tot_seq )
   {
    size_t  i;
    size_t  n;
    char   *s;
    char   *s_end;

    i = 0;
    s = seqbuff;
    s_end = seqbuff + tot_seq;
    while ( s < s_end )
       {
        if ( *s != '\0' )
            sa[i++] = s;
        s++;
       }
    n = i;
    /* printf( "Starting qsort...\n" ); */
    qsort( sa, n, sizeof( char * ), alpha_cmp );
    /* printf( "...Finished qsort\n" ); */
    return( n );
   }


size_t  common_prefix( char *a, char *b )   /* returns the length of the max */
   {                                        /* common prefix shared by a & b */
    size_t  n;

    n = 0;
    while ( *a != '\0' && *b != '\0' && *a++ == *b++ )
        n++;
    return( n );
   }


int  seq_id( char *s, char *seqbuff, char *seq2 )   /* determines which seq  */
   {                                                /* s is in, the first (0)*/
    if ( s < seqbuff )                              /* or the 2nd (1)        */
       {
        fprintf( stderr, "This should not happen: seq_id: %p\n", s );
        exit( 1 );
       }
    if ( s < seq2 )
        return( 0 );
    else
        return( 1 );
   }


size_t  seq_offset( char *s, char *seqbuff, char *seq2 )  /* returns the loc.*/
   {                                                      /* of s in its seq */
    if ( s < seqbuff )
       {
        fprintf( stderr, "This should not happen: seq_offset %p\n", s );
        exit( 1 );
       }
    if ( s < seq2 )
        return( (size_t) (s - seqbuff) );
    else
        return( (size_t) (s - seq2) );
   }

int  is_repeat( char *s, int n )
   {
    int  i, m;

    for ( i = 0; i < n_repeats; i++ )
       {
        m = strlen( repeats[i] );
        if ( n < m ) 
             m = n;
        if ( common_prefix( s, repeats[i] ) >= m )
           {
            /* printf( "@@@@found repeat\n" ); */
            return( 1 );
           }
       }
    return( 0 );
   }

int  good_seq( char *s, int n)
   {
    int    n_xs = 0;
    int    i;

    for ( i = 0; i < n; i++ )
        if ( s[i] == 'x' || s[i] == 'X' )
            n_xs++;
    return( ((double) n_xs / (double) n ) < 0.1 );
   }

void find_max_common( char    *sa[], 
                      size_t   n_sa, 
                      char    *seqbuff, 
                      char    *seq2,
                      size_t  *maxprefix,    /* output */
                      size_t  *maxind        /* output */
                     )
   {
    size_t i;
    size_t n;
    int    id, last_id;
    
    *maxprefix = 0;
    *maxind = -1;
    last_id = seq_id( sa[0], seqbuff, seq2 );
    
    for ( i = 1; i < n_sa; i++ )
       {
        id = seq_id( sa[i], seqbuff, seq2 );
        if ( id != last_id )
           {
            n = common_prefix( sa[i-1], sa[i] );
            if ( n > *maxprefix && good_seq( sa[i], n ) 
                 && ! is_repeat( sa[i], n ) )
               {
                *maxprefix = n;
                *maxind = i;
               }
           }
        last_id = id;
       }
   }



void  output_seq ( char *s, size_t n )   /* print sequence neatly */
   {
    size_t i = 0;

    while ( n-- > 0 && *s != '\0' )
       {
        putchar( *s++ );
        if ( ++i % 60 == 0 )
            putchar( '\n' );
       }
    if ( i % 60 != 0 )
        putchar( '\n' );
   }



                                /****************/
                                /* Main Program */
                                /****************/


main( int argc, char **argv )
   {
    char     file1[MAX_FILENAME+1];
    char     file2[MAX_FILENAME+1];
    char     repeatfile[MAX_FILENAME+1];
    char    *seqbuff;
    char    *seq2;
    size_t   tot_seq;
    char   **sa;
    size_t   n_sa;
    size_t   max_len_f, max_len_r;
    size_t   max_ind_f, max_ind_r;
    size_t   off1, off2;
    size_t   max_len;
    char     dir;
    size_t   len1, len2;

    parse_args( argc, argv, file1, file2, repeatfile );
    /* printf( "File 1 is %s  %d\n", file1, filesize( file1 ) );
       printf( "File 2 is %s  %d\n", file2, filesize( file2 ) ); */

    if ( *repeatfile != '\0' )
        load_repeats( repeatfile );

    tot_seq = filesize( file1 ) + filesize( file2 ) + 2;
    if ( !( seqbuff = (char *) malloc( tot_seq )) )
       {
        fprintf( stderr, "unable to malloc %d bytes for sequences\n", tot_seq);
        exit( 1 );
       }

    len1 = load_seq( file1, seqbuff, tot_seq );
    seq2 = seqbuff + len1 + 1;
    len2 = load_seq( file2, seq2, tot_seq - (len1+1));
    tot_seq = len1 + len2;    
    /* printf( "tot_seq is %d\n", tot_seq ); */

    if ( !( sa = (char **) malloc( tot_seq * sizeof(char *)  )) )
       {
        fprintf( stderr, "unable to malloc %d bytes for suffix array\n", 
                         tot_seq * sizeof(char *) );
        exit( 1 );
       }
 
    max_len = 0;
    n_sa = create_suffix_array( sa, seqbuff, tot_seq );
    find_max_common( sa, n_sa, seqbuff, seq2, &max_len_f, &max_ind_f );
    /* printf( "foward max %d ind %d\n", max_len_f, max_ind_f ); */
    if ( max_len_f > 0 )
       {
        max_len = max_len_f;
        dir = 'f';
        if ( seq_id( sa[max_ind_f], seqbuff, seq2 ) == 1 )
           {
            off1 = seq_offset( sa[max_ind_f - 1], seqbuff, seq2 );
            off2 = seq_offset( sa[max_ind_f], seqbuff, seq2 );
           }
        else if ( seq_id( sa[max_ind_f - 1], seqbuff, seq2 ) == 1 )
           {
            off2 = seq_offset( sa[max_ind_f - 1], seqbuff, seq2 );
            off1 = seq_offset( sa[max_ind_f], seqbuff, seq2 );
           }
        else
           {
            fprintf( stderr, "this should not happen: bad seq_id f\n" );
            exit( 1 );
           }
       }
/*    else
        printf( "No forward overlap\n" ); */

    /* Now: reverse-complement the 2nd sequence and repeat the whole */
    /* procedure                                                     */

    rc_in_situ( seq2 );

    n_sa = create_suffix_array( sa, seqbuff, tot_seq );
    find_max_common( sa, n_sa, seqbuff, seq2, &max_len_r, &max_ind_r );

    if ( max_len_r > max_len )
       {
        max_len = max_len_r;
        dir = 'r';
        if ( seq_id( sa[max_ind_r], seqbuff, seq2 ) == 1 )
           {
            off1 = seq_offset( sa[max_ind_r - 1], seqbuff, seq2 );
            off2 = seq_offset( sa[max_ind_r], seqbuff, seq2 );
           }
        else if ( seq_id( sa[max_ind_r - 1], seqbuff, seq2 ) == 1 )
           {
            off2 = seq_offset( sa[max_ind_r - 1], seqbuff, seq2 );
            off1 = seq_offset( sa[max_ind_r], seqbuff, seq2 );
           }
        else
           {
            fprintf( stderr, "this should not happen: bad seq_id f\n" );
            exit( 1 );
           }
       }
/*    else
        printf( "No forward overlap\n" ); */

    if ( max_len > 0 )
       {
        if ( print_seq )
           printf( ">" );
        printf( "%8d %6.2f %6.2f %8d %c %8d %8d %8d\n", max_len, 
                                     (100.0 * max_len ) / len1,
                                     (100.0 * max_len ) / len2,
                                     off1, 
                                     dir, off2, len1, len2 );
        if ( print_seq )
            output_seq( seqbuff + off1, max_len );
       }
    else
       printf( "no overlap\n" ); 

    free( seqbuff );
    free( sa );
    exit( 0 );
   }



