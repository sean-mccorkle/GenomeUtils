/* File:         suffixtree.c                                                */
/* Programmer:   Sean R. McCorkle                                            */
/* Language:     C                                                           */
/* Description:  complete library for manipulating suffix trees for DNA and  */
/*               amino acid sequences                                        */
/*                                                                           */
/* $Id$  */
/*****************************************************************************/

#include <stdio.h>
#include <string.h>
#include "suffixtree.h"
#include "io.h"
                        /**************************/
                        /* Suffix Tree Operations */
                        /**************************/


ST_NODE *new_st_node( char c, int loc, ST_NODE *children, ST_NODE *next )
   {
    ST_NODE *s;

    s = (ST_NODE *) malloc_safely( sizeof( ST_NODE ) );
    s->c = c;
    s->loc = loc;
    s->children = children;
    s->next = next;
    return( s );
   }


ST_NODE *linear_suffix_tree( char *s, int loc )
   {
    if ( *s == '\0' )
        return( NULL );
    else
        return( new_st_node( *s, 
                             (*s == '$' ? loc : 0 ), 
                             linear_suffix_tree( s+1, loc ),
                             NULL
                            )
              );
   }


ST_NODE *search_st_list( ST_NODE *l, char c )
   {
    while ( l != NULL && l->c != c )
        l = l->next;
    return( l );
   }

/*    add_to_suffix_tree( t, s ); */
/* adds substring s to existing tree */

void  add_to_suffix_tree( ST_NODE *t, char *s, int loc )
   {
    ST_NODE *p;

    if ( s[0] == '\0' )
        return;
    if ( p = search_st_list( t->children, s[0] ) )/* if we find a lower node */
        add_to_suffix_tree( p, s+1, loc );        /* follow & recurse        */
    else                            
       { 
        p = linear_suffix_tree( s, loc );       /* if not, then we add a new */
        p->next = t->children;                  /* linear sub-branch to the  */
        t->children = p;                        /* children of the curr. node*/
       }
    return;
   }


/* tree = create_suffix_tree( s ); */

ST_NODE *create_suffix_tree( char *s )
   {
    ST_NODE *t;
    int      i = 0;

    t = new_st_node( '^', 0, NULL, NULL );
    while ( *s != '\0' && *s != '$' )
        add_to_suffix_tree( t, s++, ++i );
    return( t );
   }


void  print_suffix_tree( ST_NODE *t, int indent )
   {
    int  i;
    ST_NODE *c;

    if ( t == NULL )
        return;
    for ( i = 0; i < indent; i++ )
        putchar( ' ' );
    printf( "%c", t->c );
    if ( t->loc > 0 )
        printf( "%d", t->loc );
    printf( "\n" );

    print_suffix_tree( t->children, indent + 4 );
    print_suffix_tree( t->next, indent );
   }


void  enumerate_st_hits( ST_NODE *t )
   {
    if ( t == NULL )
        return;
    if ( t->c == '$' )
        printf( "%d\n", t->loc );
    else    
        enumerate_st_hits( t->children );
    enumerate_st_hits( t->next );
   }



void  search_suffix_tree( ST_NODE *t, char *s )
   {
    ST_NODE *p;

    if ( s[0] == '\0' )
       {
        enumerate_st_hits( t->children );
        return;
       }
    if ( !( p = search_st_list( t->children, s[0] ) ) )
        return;                                         /* not successful */
    search_suffix_tree( p, s+1 );
   }





