
#include <stdio.h>

void test_routine( char *a, char *b );

main()
   {
    printf( "hello\n" );
    test_routine( "x" );
   }

void test_routine( char *a, char *b )
   {
    if ( a )
       printf( "a [%s]\n", a );
    if ( b )
       printf( "b [%s]\n", b );
   }
