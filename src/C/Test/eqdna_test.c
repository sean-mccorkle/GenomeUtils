#include <stdio.h>
#include "../seqlib.h"


char lets[] = { 'A',   'C',   'G',   'T',
                'M',   'R',   'W',   'S',
                'Y',   'K',   'V',   'H',
                'D',   'B',   'N',   'X', '\0' };

main()
   {
    int  i, j;

    for ( i = 0; lets[i]; i++ )
        for ( j = 0; lets[j]; j++ )
            printf( "%c  %c  %s\n", lets[i], lets[j], 
                                    eqdna( lets[i], lets[j] ) ? "yes" : "no" 
                  ); 

   }
