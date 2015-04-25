
#include <stdio.h>
#include <stdlib.h>
#include "io.h"

void  bailout( char *s )
   {
    fprintf( stderr, "%s\n", s );
    exit( 2 );
   }

void  *malloc_safely( size_t size )
   {
    void *p;

    if ( ! (p = malloc( size ) ) )
        bailout( "failed to malloc new st_node" );
    return( p );
   }






