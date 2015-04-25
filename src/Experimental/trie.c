/* File:         trie.c                                                      */
/* Programmer:   Sean R. McCorkle                                            */
/* Language:     C                                                           */
/* Description:  library for manipulating tries for DNA and amino acid       */
/* sequences                                                                 */
/*                                                                           */
/* $Id$  */
/*****************************************************************************/

#include <assert.h>
#include <stdio.h>
#include <string.h>
#include "trie.h"
#include "io.h"

#define DFS_STACK_SIZE 10000    /* stack limit for DFSprocess() */

void      add_trie_child( TR_NODE *t, TR_NODE *ch );
void      add_trie_leaf( TR_NODE *t, TR_LEAF *l );
void      print_leaves( TR_LEAF *l, int indent, void (*print_data)(void *) );
TR_LEAF  *new_tr_leaf( void *data, TR_LEAF *next );
TR_NODE  *new_tr_node( char c, TR_LEAF *l, TR_NODE *children, TR_NODE *next );
TR_NODE  *search_tr_list( TR_NODE *l, char c );
TR_NODE  *trie_limb( char *s, void *data );
char     *DFSpath( void );

typedef enum { FORWARD = 'F', REVERSE = 'R' } DIRECTION;
 
typedef struct {
                 int             string_id;
                 int             pos;
                 DIRECTION       dir;
               } POS_REC;


POS_REC  *pos_rec( int id, int pos, DIRECTION d )
   {
    POS_REC *p;

    p = (POS_REC *) malloc_safely( sizeof( POS_REC ) );
    p->string_id = id;
    p->pos = pos;
    p->dir = d;
    return( p );
   }

void print_pos( void *d )
   {
    POS_REC *p = (POS_REC *) d;
    printf( "  (%d %d %c)\n", p->string_id, p->pos, p->dir );
   }

void analyze( TR_LEAF *l )
   {
    if ( l == NULL )
       return;
    printf( "%s\n", DFSpath() );
    print_leaves( l, 5, print_pos );
   }

main()
   {
    TR_NODE *t;
 
    t = new_tr_node( '#', NULL, NULL, NULL );
    add_trie_string( t, "acccgt", pos_rec( 1, 25, FORWARD ) );
    add_trie_string( t, "gta", pos_rec( 2, 50, REVERSE ) );
    add_trie_string( t, "accgt", pos_rec( 3, 75, FORWARD ) );
    add_trie_string( t, "accgt", pos_rec( 4, 100, REVERSE ) );
    add_trie_string( t, "gta", pos_rec( 5, 125, FORWARD ) );
    print_trie( t, 0, print_pos );
    DFSprocess( t, analyze );
   }


/**************************************************/
/* procedures & functions intended for public use */
/**************************************************/

/* add in string s with first character starting as children of t */


void add_trie_string( TR_NODE *t, char *s, void *data )
   {
    TR_NODE *ch;

    if ( *s == '\0' ) /* string is exactly in trie.  just add leaf */
        add_trie_leaf( t, new_tr_leaf( data, NULL ) );
    else if ( (ch = search_tr_list( t->children, *s )) != NULL )
        add_trie_string( ch, s+1, data );   /* recurse */
    else
        add_trie_child( t, trie_limb( s, data ) );
   }


/* process all branches, Depth First Search */

void  DFSprocess( TR_NODE *t, void (*proc)(TR_LEAF *) )
   {
    if ( t == NULL )
        return;
    DFSpush( t->c );
    if ( proc != NULL )
        (*proc)( t->leaves );
    DFSprocess( t->children, proc );
    DFSpop();
    DFSprocess( t->next, proc );
   }

/* print out suffix tree */

void  print_trie( TR_NODE *t, int indent, void (* print_data)(void *) )
   {
    int  i;

    if ( t == NULL )
        return;
    for ( i = 0; i < indent; i++ )
        putchar( ' ' );
    putchar( t->c );
    printf( "\n" );
    print_leaves( t->leaves, indent, print_data );

    print_trie( t->children, indent + 4, print_data );
    print_trie( t->next, indent, print_data );
   }

void  print_leaves( TR_LEAF *l, int indent, void (* print_data)(void *) )
   {
    int      i;
    if ( l == NULL )
        return;
    for ( i = 0; i < indent; i++ )
        putchar( ' ' );
    if ( print_data != NULL )
        (*print_data)( l->data );
    print_leaves( l->next, indent, print_data );
   }

/********************************************************/
/* procedures & functions intended only for private use */
/********************************************************/


TR_NODE *trie_limb( char *s, void *data )
   {
    assert( *s != '\0' ); /* add_trie_string should never allow this to be 0 */
    if ( *(s+1) != '\0' )   /* recurse unless string is one char long */
        return( new_tr_node( *s, NULL, trie_limb( s+1, data ), NULL ) );
    else                   /* only one char left in string. */
        return( new_tr_node( *s, new_tr_leaf( data, NULL ), NULL, NULL ) );
   }

/* adds c subtree to t's child list */

void add_trie_child( TR_NODE *t, TR_NODE *ch )
   {
    assert( t != NULL && ch != NULL && ch->next == NULL );
    ch->next = t->children;
    t->children = ch;
   }

void  add_trie_leaf( TR_NODE *t, TR_LEAF *l )
   {
    assert( t != NULL && l != NULL && l->next == NULL );
    l->next = t->leaves;
    t->leaves = l;
   }


TR_NODE *search_tr_list( TR_NODE *l, char c )
   {
    while ( l != NULL && l->c != c )
        l = l->next;
    return( l );
   }

/* globals for use only by DFSpush(), DFSpop(), DFSpath() */

static char dfs_stack[DFS_STACK_SIZE];
static int  dfs_sp = 0;

void  DFSpush( char c ) 
   {
    if ( dfs_sp >= DFS_STACK_SIZE )
       {
        fprintf( stderr, "DFS stack overflow\n" );
        exit( 1 );
       }
    dfs_stack[dfs_sp++] = c;
    dfs_stack[dfs_sp] = '\0';
   }

void  DFSpop( void ) 
   {
    assert( dfs_sp > 0 );
    dfs_stack[--dfs_sp] = '\0';
   }

char *DFSpath( void )
   { 
     return( dfs_stack + 1 );   /* root node is '#' - don't forget! */
   }


TR_NODE *new_tr_node( char c, TR_LEAF *leaf, TR_NODE *children, TR_NODE *next )
   {
    TR_NODE *n;

    n = (TR_NODE *) malloc_safely( sizeof( TR_NODE ) );
    n->c = c;
    n->leaves = leaf;
    n->children = children;
    n->next = next;
    return( n );
   }

TR_LEAF  *new_tr_leaf( void *data, TR_LEAF *next )
   {
    TR_LEAF *l;

    l = (TR_LEAF *) malloc_safely( sizeof( TR_LEAF ) );
    l->data = data;
    l->next = next;
    return( l );
   }



