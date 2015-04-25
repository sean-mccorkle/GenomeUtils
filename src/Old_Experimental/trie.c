/* File:         trie.c                                                      */
/* Programmer:   Sean R. McCorkle                                            */
/* Language:     C                                                           */
/* Description:  library for manipulating tries for DNA and amino acid       */
/* sequences                                                                 */
/*                                                                           */
/* $Id$  */
/*****************************************************************************/

#include <stdio.h>
#include <string.h>
#include "trie.h"
#include "io.h"

void      add_trie_child( TR_NODE *t, TR_NODE *ch );
void      add_trie_leaf( TR_NODE *t, TR_LEAF *l );
TR_LEAF  *new_tr_leaf( int id, int pos, DIRECTION d, TR_LEAF *next );
TR_NODE  *new_tr_node( char c, TR_LEAF *l, TR_NODE *children, TR_NODE *next );
TR_NODE  *search_tr_list( TR_NODE *l, char c );
TR_NODE  *trie_limb( char *s, int string_id, int pos, DIRECTION d );


main()
   {
    TR_NODE *t;
 
    t = new_tr_node( '#', NULL, NULL, NULL );
    add_trie_string( t, "acccgt", 1, 25, FORWARD );
    add_trie_string( t, "gta", 2, 50, REVERSE  );
    add_trie_string( t, "accgt", 3, 75, FORWARD );
    printf( "next two should already be in\n" );
    add_trie_string( t, "accgt", 4, 100, REVERSE );
    add_trie_string( t, "gta", 5, 125, FORWARD );
    print_trie( t, 0 );
   }


/**************************************************/
/* procedures & functions intended for public use */
/**************************************************/

/* add in string s with first character starting as children of t */

void add_trie_string( TR_NODE *t, char *s, int string_id, int pos, 
                      DIRECTION d )
   {
    TR_NODE *ch;

    if ( *s == '\0' )
       {
        /* string is exactly in trie.  add stuff */
        printf( "this string already in at %c\n", t->c );
        return;
       }
    if ( (ch = search_tr_list( t->children, *s )) != NULL )
        add_trie_string( ch, s+1, string_id, pos, d );   /* recurse */
    else
        add_trie_child( t, trie_limb( s, string_id, pos, d ) );
   }

/* print out suffix tree */

void  print_trie( TR_NODE *t, int indent )
   {
    int  i;

    if ( t == NULL )
        return;
    for ( i = 0; i < indent; i++ )
        putchar( ' ' );
    putchar( t->c );
    printf( "\n" );

    print_trie( t->children, indent + 4 );
    print_trie( t->next, indent );
   }

/********************************************************/
/* procedures & functions intended only for private use */
/********************************************************/


TR_NODE *trie_limb( char *s, int string_id, int pos, DIRECTION d )
   {
    if ( *s == '\0' )  /* add_trie_string() should never invoke this with \0 */
       {
        fprintf( stderr, "error: null string passed to trie_limb" );
        exit( 1 );
       }
    if ( *(s+1) != '\0' )   /* recurse unless string is one char long */
        return( new_tr_node( *s, NULL,
                             trie_limb( s+1, string_id, pos, d ), NULL 
                           ) 
              );
    else     /* only one char left in string. */
        return( new_tr_node( *s, make_tr_leaf( string_id, pos, d, NULL ),
                              NULL, NULL 
                            )
              );
   }

/* adds c subtree to t's child list */

void add_trie_child( TR_NODE *t, TR_NODE *ch )
   {
    if ( ch->next != NULL )
       {
        fprintf(stderr,"error - got non-null next in add_child: (%c)\n",ch->c);
        exit( 0 );
       }
    ch->next = t->children;
    t->children = ch;
   }

void  add_trie_leaf( TR_NODE *t, TR_LEAF *l )
   {
    if ( t == NULL )
       {
        fprintf( stderr, "error - add_trie_leaf got null t\n" );
        exit( 0 );
       }
    if ( l->next != NULL )
       {
        fprintf(stderr,"error - got non-null next in add_trie_leaf\n" );
        exit( 0 );
       }
    l->next = t->
   }


TR_NODE *search_tr_list( TR_NODE *l, char c )
   {
    while ( l != NULL && l->c != c )
        l = l->next;
    return( l );
   }


TR_NODE *new_tr_node( char c, TR_LEAF *leaf, TR_NODE *children, TR_NODE *next )
   {
    TR_NODE *n;

    n = (TR_NODE *) malloc_safely( sizeof( TR_NODE ) );
    n->c = c;
    n->leaf = leaf;
    n->children = children;
    n->next = next;
    return( n );
   }

TR_LEAF  *new_tr_leaf( int id, int pos, DIRECTION d, TR_LEAF *next )
   {
    TR_LEAF *l;

    l = (TR_LEAF *) malloc_safely( sizeof( TR_LEAF ) );
    l->string_id = id;
    l->pos = pos;
    l->dir = d;
    l->next = next;
    return( l );
   }


