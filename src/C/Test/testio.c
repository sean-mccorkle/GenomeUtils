
#include <stdio.h>
#include "seqlib.h"

char *progname;

main( int argc, char **argv )
   {
    int         s;
    static char hdr[256];
    static char seq[MAX_SEQ_LENGTH];

    while ( --argc > 0 )
       {
        argv++;
        if ( (s = open_seqf( *argv, STD_ALLOWED, STD_IGNORED ) ) < 0 )
            fprintf( stderr, "can't open %s\n", *argv );
        else
           {
            while ( read_seq( s, 256, hdr, MAX_SEQ_LENGTH, seq ) >= 0 )
               {
                printf( ">%s\n", hdr );
                print_seq( stdout, seq, 0, 0, 0, 0, 0, 0 );
               }
            close_seqf( s );
           }
       }
   }

