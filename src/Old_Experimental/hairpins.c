#include <stdio.h>
#include <string.h>
#include "io.h"
#include "suffixtree.h"

#define  SMAX  100000
#define  STACKMAX 1000

void  getstr( char *s, int size )
   {
    int l;

    if ( ! fgets( s, size, stdin ) )
        exit( 2 );
    s[ strlen( s ) - 1 ] = '$';   /* blow away carraige return */
   }

void reverse_complement( char *r, char *s )
   {
    strlen
   }

main( int argc, char **argv )

   {
    static char s[SMAX], r[SMAX];
    ST_NODE     *f_tree, *r_tree;
    static char st[STACKMAX];
    int         depth;
    int         len;

    getstr( s, SMAX );
    printf( "hairpins: [%s]\n", s );    
    f_tree = create_suffix_tree( s );
    printf( "suffix tree for [%s]\n", s );
    print_suffix_tree( f_tree, 0 );

   }

