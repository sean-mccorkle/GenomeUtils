#include <errno.h>
#include <stdio.h>
#include <string.h>
#include "seqlib.h"
#include "trie.h"

#define  MAX_SEQ_LEN  50000
#define  MAX_FILE_LEN 256

char filename[MAX_FILE_LEN];

typedef enum { FORWARD = 'F', REVERSE = 'R' } DIRECTION;
 
typedef struct {
                 int             pos;
                 DIRECTION       dir;
               } POS_REC;

POS_REC  *pos_rec( int pos, DIRECTION d )
   {
    POS_REC *p;

    p = (POS_REC *) malloc_safely( sizeof( POS_REC ), "pos_rec" );
    p->pos = pos;
    p->dir = d;
    return( p );
   }


void  print_pos( void *d )
   {
    POS_REC *p = (POS_REC *) d;
    printf( "  (%d %c)\n", p->pos, p->dir );
   }


void  analyze( TR_LEAF *l )
   {
    int      n_for;
    int      n_rev;
    TR_LEAF *q;

    if ( l == NULL )
        return;
    for (  q = l, n_for = n_rev = 0;  q != NULL;  q = q->next  )
       {
        if ( ((POS_REC *)q->data)->dir == FORWARD )
            n_for++;
        else
            n_rev++;
       }

    if ( n_for + n_rev > 1 )
       {
        printf( "%s: %d %d\n", DFSpath(), n_for, n_rev );
        for ( q = l; q != NULL; q = q->next )
           print_pos( q->data );
        printf( "----\n" );
       }
   }


void  parse_args( int argc, char **argv )
   {
    ++argv; --argc;
    if ( argc != 1 )
       {
        fprintf( stderr, "need a file name\n" );
        exit( 1 );
       }
    strncpy( filename, *argv, MAX_FILE_LEN );
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

void  DFSsearchmax( TR_NODE *t )
   {
    static int  max_depth = 0;

    if ( t == NULL )
        return;
    if ( t->count < 2 )
        return;
    DFSpush( t->c );
    DFSsearchmax( t->children );
    if ( DFSdepth() >= max_depth )
       {
        max_depth = DFSdepth();
        printf( "%d: %s\n", max_depth, DFSpath() );
       }
    DFSpop();
    DFSsearchmax( t->next );
   }


main( int argc, char **argv )
   {
    int            seqf;
    static char    seq[MAX_SEQ_LEN];
    static char    rev[MAX_SEQ_LEN];
    int            len;  /* length of seq */
    TR_NODE       *t;
    char          *p;

    parse_args( argc, argv );
    printf( "File is %s\n", filename );

    if ( ( seqf = open_seqf( filename, STD_ALLOWED, STD_IGNORED ) < 0 ) )
       {
        fprintf( stderr, "%s: fatal; ", progname );
        perror( filename );
        exit( errno );
       }
    if ( get_seq( seqf, seq, &len ) <= 0 )
       {
        fprintf( stderr, "%s: no sequence in file %s\n", progname, filename );
        exit( 1 );
       }
    printf( "seq len is %d\n", len );

    dna_reverse_comp( seq, rev );

    printf( "Building suffix tree...\n" );
    t = create_trie();
    printf( "...forward\n" );
    for ( p = seq; *p != '\0'; p++ )
        add_trie_string( t, p, pos_rec( (int)( p - seq ), FORWARD ) );

    printf( "...reverse\n" );
    for ( p = rev; *p != '\0'; p++ )
        add_trie_string( t, p, pos_rec( (int)( p - rev ), REVERSE ) );

#ifdef NO
    print_trie( t, 0, print_pos );  
#endif
    printf( "...Done.\nBeginning DFS analysis of tree....\n" );

    DFScount( t );
    /* print_trie( t, 0, print_pos );   */
    /* DFSprocess( t, analyze ); */
    DFSsearchmax( t );
   }




