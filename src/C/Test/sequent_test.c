#include <stdio.h>
#include "../seqlib.h"


main()
   {
    char c, C;
    for ( c = 'a', C = 'A'; c <= 'z'; c++, C++ )
        printf( "%c: %2d     %c: %2d\n", c, SEQUENT(c), C, SEQUENT(C) );
   }
