#include <errno.h>
#include <stdio.h>
#include <string.h>
#include "seqlib.h"
#include "trie.h"

#define  MAX_SEQ_LEN  50000
#define  MAX_FILE_LEN 256

int  min_length = 12;

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

void  print_pos2( void *d )
   {
    POS_REC *p = (POS_REC *) d;
    printf( "  (%d %c)", p->pos, p->dir );
   }

void  print_pos3( void *d )
   {
    POS_REC *p = (POS_REC *) d;
    printf( " %c%d",  p->dir, p->pos);
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
    int            ret;
    static char    hdr[MAX_HDR_LENGTH];

    ret = read_seq( seqf, MAX_HDR_LENGTH, hdr, MAX_SEQ_LENGTH, seq );
    if ( ret > 0 ) 
       {
        *len = strlen( seq );
        uppercase( seq, seq );
       }
    return( ret );
   }

void  output_seq( TR_NODE *t )
   {
    static int  number = 0;
    char       *s;
    int         i, n;

    printf( ">rep%d %d", ++number, DFSdepth() );
    DFSprintleaves( t, print_pos3 );
    printf( "\n" );

    s = DFSpath();
    n = strlen( s );
    for ( i =  60; i <= n; i += 60, s += 60 )
        printf( "%60.60s\n", s );
    if ( *s != '\0' )
        printf( "%s\n", s );
   }

void  DFSrepeats( TR_NODE *t )
   {
    if ( t == NULL )
        return;
    if ( t->count < 2  )
       {
        if ( DFSdepth() >= min_length )
            output_seq( t );
        return;
       }

    DFSpush( t->c );
    if ( DFSdepth() >= min_length && t->leaves != NULL )
        output_seq( t );
    DFSrepeats( t->children );
    DFSpop();
    DFSrepeats( t->next );
   }


main( int argc, char **argv )
   {
    int            seqf;
    static char    seq[MAX_SEQ_LEN];
    static char    rev[MAX_SEQ_LEN];
    static char    tmp[MAX_SEQ_LEN];
    int            len;  /* length of seq */
    TR_NODE       *t;
    char          *p;
    int            i;

    int            trunc = 500;

    parse_args( argc, argv );
    fprintf( stderr, "File is %s\n", filename );

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
    fprintf( stderr, "seq len is %d\n", len );

    dna_reverse_comp( seq, rev );

    fprintf( stderr, "Building suffix tree...\n" );
    t = create_trie();
    fprintf( stderr, "...forward\n" );
    for ( p = seq, i = 1; *p != '\0'; p++ )
       {
        strncpy( tmp, p, trunc );
        add_trie_string( t, tmp, pos_rec( (int)( p - seq ), FORWARD ) );
        if ( (++i % 1000) == 0 )
            fprintf( stderr, "...%d\n", i );
       }

    fprintf( stderr, "...reverse\n" );
    for ( p = rev, i = 0; *p != '\0'; p++ )
       {
        strncpy( tmp, p, trunc );
        add_trie_string( t, tmp, pos_rec( (int)( p - rev ), REVERSE ) );
        if ( (++i % 1000) == 0 )
            fprintf( stderr, "...%d\n", i );
       }

#ifdef NO
    print_trie( t, 0, print_pos );  
#endif
    fprintf( stderr, "...Done.\nBeginning DFS analysis of tree....\n" );

    DFScount( t );
    /* print_trie( t, 0, print_pos );   */
    /* DFSprocess( t, analyze ); */
    DFSrepeats( t );
    fprintf( stderr, "...Done.\nProgram ends.\n" );
   }




