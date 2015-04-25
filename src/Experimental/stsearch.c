#include <stdio.h>
#include <string.h>
#include "io.h"
#include "suffixtree.h"

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
    printf( "[%s]\n", s );    
    tree = create_suffix_tree( s );
    printf( "suffix tree for [%s]\n", s );
    print_suffix_tree( tree, 0 );
    gets( a );
    printf( "searching for [%s]\n", a );
    search_suffix_tree( tree, a );
   }






