#include <stdio.h>
#include <string.h>

struct st_node {
                char   c;
                int    loc;
                struct st_node *children;
                struct st_node *next;
               };
typedef struct st_node  ST_NODE;


void  bail( char *s )
   {
    fprintf( stderr, "%s\n", s );
    exit( 2 );
   }

                        /**************************/
                        /* Suffix Tree Operations */
                        /**************************/


ST_NODE *new_st_node( void )
   {
    ST_NODE *s;

    if ( ! (s = (ST_NODE *) malloc( sizeof( ST_NODE ) ) ) )
        bail( "failed to malloc new st_node" );
    s->c = ' ';
    s->loc = 0;
    s->children = NULL;
    s->next = NULL;
    return( s );
   }


ST_NODE *linear_suffix_tree( char *s, int loc )
   {
    ST_NODE *t;

    if ( *s == '\0' )
        return( NULL );
    t = new_st_node();
    t->c = s[0];
    if ( *s == '$' )
        t->loc = loc;
    t->children = linear_suffix_tree( s+1, loc );
    return( t );
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
        add_to_suffix_tree( p, s+1, loc );      /* follow & recurse        */
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

    t = new_st_node();
    t->c = '^';
    while ( *s != '\0' )
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
        enumerate_st_hits( t );
        return;
       }
    if ( !( p = search_st_list( t->children, s[0] ) ) )
        return;                                         /* not successful */
    search_suffix_tree( p, s+1 );
   }

#define  SMAX  100

void getstr( char *s, int size )
   {
    int l;

    if ( ! fgets( s, size, stdin ) )
        exit( 2 );
    s[ strlen( s ) - 1 ] = '$';   /* blow away carraige return */
   }


main( int argc, char **argv )

   {
    static char s[SMAX], a[SMAX];
    ST_NODE     *tree;

    getstr( s, SMAX );
    
    tree = create_suffix_tree( s );
    /*printf( "suffix tree for [%s]\n", s );*/
    /*print_suffix_tree( tree, 0 ); */
    gets( a );
    printf( "searching for [%s]\n", a );
    search_suffix_tree( tree, a );
   }






