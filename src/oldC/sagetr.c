#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifndef SOLARIS        /* Sun's are really a pain in the neck sometimes */
#include <unistd.h>
#endif
#include "seqlib.h"
#include "trie.h"

#define  MAX_SEQ_LEN  200000

typedef enum { FORWARD = 'F', REVERSE = 'R' } DIRECTION;
typedef enum { TRUE = 1, FALSE = 0 } BOOLEAN;
 
typedef struct {
                 int             string_id;
                 int             pos;
                 DIRECTION       dir;
                 BOOLEAN         real_tag;  /* F => not the real tag */
               } POS_REC;

/* globals which are potentially set by command line options */

char  *re_seq     = "CATG";
int    tag_extent = 10;
int    n_files;
char **files;

/* globals for calculations */

int    n_tags = 0;
int    n_mult_tags = 0;

POS_REC  *pos_rec( int id, int pos, DIRECTION d, BOOLEAN real_tag )
   {
    POS_REC *p;

    p = (POS_REC *) malloc_safely( sizeof( POS_REC ), "pos_rec" );
    p->string_id = id;
    p->pos = pos;
    p->dir = d;
    p->real_tag = real_tag;
    return( p );
   }


void  print_pos( void *d )
   {
    POS_REC *p = (POS_REC *) d;
    printf( "  (%d %d %c %d)\n", p->string_id, p->pos, p->dir, p->real_tag );
   }


void  analyze( TR_LEAF *l )
   {
    int      n;
    int      n_real;
    TR_LEAF *q;

    if ( l == NULL )
        return;
    for (  q = l, n = 1, n_real = 0 ;  q != NULL;  q = q->next, n++  )
       {
        if ( ((POS_REC *)(q->data))->real_tag )
            n_real++;
       }

    if ( n_real > 0 )
        n_tags++;

    if ( n_real > 1 )
       {
        n_mult_tags++;
        printf( "%s: %d/%d\n", DFSpath(), n_real, n );
#ifdef NO
        for ( q = l; q != NULL; q = q->next )
           print_pos( q->data );
        printf( "----\n" );
#endif
       }
   }


void  parse_args( int argc, char **argv )
   {
    int          c;
    int          i;
    extern char *optarg;
    extern int   optind;
          
    while ( (c = getopt(argc, argv, "e:s:h") ) != -1 )
        switch (c) 
           {
            case 'e':
                       tag_extent = atoi( optarg );
                       break;
            case 's':
                       re_seq = strdup( optarg );
                       uppercase( re_seq, re_seq );
                       break;
            case 'h': 
                       printf( "no help yet. sorry.\n" );
                       exit( 0 ); 
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


int  get_seq( int seqf, char *seq, int *len )
   {
    int            n_hdr_lines;
    char         **hdr;
    int            ret;

    load_header( seqf, &n_hdr_lines, &hdr );
    free_lines( n_hdr_lines, hdr );
    ret = load_seq( seqf, MAX_SEQ_LEN, seq, len );
    if ( ret > 0 ) 
        uppercase( seq, seq );
    return( ret );
   }


int  main( int argc, char **argv )
   {
    TR_NODE       *t;
    int            seqf;
    int            i;
    static char    seq[MAX_SEQ_LEN];
    int            len;  /* length of seq */
    int            id;   /* seq id (incrementing counter) */
    char          *p;    /* position of restriction site in seq */
    static char    tag[200];
    int            tag_len;
    int            negative_fs = 0;  /* count of seqs with no re site forward*/
    POS_REC       *r;

    parse_args( argc, argv );
    printf( "tag extent: %d\n", tag_extent );
    printf( "RE seq: %s\n", re_seq );
    tag_len = strlen( re_seq ) + tag_extent;
    printf( "tag len: %d\n", tag_len );
    t = create_trie();
    id = 0;
    printf( "files:\n" );
    for ( i = 0; i < n_files; i++ )
       {
        printf( "   %s\n", files[i] );
        if ( ( seqf = open_seqf( files[i], STD_ALLOWED, STD_IGNORED ) < 0 ) )
           {
            fprintf( stderr, "%s: fatal; ", progname );
            perror( files[i] );
            exit( errno );
           }
        while ( get_seq( seqf, seq, &len ) > 0 )
           {
            strcat( seq, "AAAAAAAAAAAAAAAAAAAA" ); /* fix this! */
#ifdef DEBUG
            printf( "sequence loaded\n" );
            print_seq( stdout, seq, 0, 0, 5, 0, 0, 0 );
#endif
            ++id;
            r = NULL;
            if ( p = strstr( seq, re_seq ) )
               {
                do {
                    strncpy( tag, p, tag_len );
                    /*printf( "id %d, tag is [%s]\n", id, tag); */
                    r = pos_rec( id, (int)(p-seq), FORWARD, FALSE );
                    add_trie_string( t, tag, r );
                   }
                while ( p = strstr( p + 1, re_seq ) );
                r->real_tag = TRUE;
               }
            else
                negative_fs++;
           }
        close_seqf( seqf );
       }
#ifdef NO
    print_trie( t, 0, print_pos );  
#endif
    printf( "finished building trie\n" );
    DFSprocess( t, analyze );
    printf( "normal exit\n" );
    printf( "%d sequences.  %d negatives (forward)\n", id, negative_fs );
    printf( "%d tags, %d are multiple tags\n", n_tags, n_mult_tags );
   }





