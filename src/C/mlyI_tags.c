/* prints out MlyI tags                   */
/*          site   +5       tag           */
/*          GAGTCNNNNN] [NNNNNNNNNNNNNNNN */

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#ifndef SOLARIS        /* Sun's are really a pain in the neck sometimes */
#include <unistd.h>    /* for getopt(), BSD, Irix, Linux */
#endif


#define  MAXTAG      50
#define  MAX_HDR_LEN 256
#define  MLYI_SEQ    "GAGTC"

/***********/
/* Globals */
/***********/

int  seq_num = 0;
int  pos = 0; 

/*  GAGTCNNNNN] [NNNNNNNNNNNNNNNN */
int  tag1_len = 10;
int  tag2_len = 16;
int  mlyI_pos = 0;   /* (from 0) */
int  tag_len;       /*   = tag1_len + tag2_len; */

int  verbose = 0;      /* set by -v */
int  print_hdrs = 0;

char  hdr[MAX_HDR_LEN+1];

void  help( void )
   {
    printf( "no help.  sorry.\n" );
   }
/**************************************************************************/
/* parse_ags() - read command line arguments and options and set globals. */
/* returns remaining arguments as nfiles and filenames array              */
/* note to self: the invocations of atoi() and strdup() could probably be */
/* made safer with more checks                                            */
/**************************************************************************/

/* note to self: put in length checks on -a, -l, etc */

void  parse_args( int argc, char **argv, int *nfiles, char ***filenames )
   {
    extern char *optarg;
    extern int   optind;
    int          c;
    static char *def_files[] = { "-", "" };

    while ( (c = getopt( argc, argv, "v" ) ) != -1 )
        switch( c )
           {
            case  'v':   verbose = 1;
                         break;
            default:     help();
                         exit( 1 );
           }
    argc -= optind;
    argv += optind;
    if ( argc > 0 )
       {
        *nfiles = argc;
        *filenames = argv;
       }
    else
       {
        *nfiles = 1;
        *filenames = def_files;
       }
   }

/*********************************************************************/
/* open_file() opens a file or returns stdin if name is "-", or does */
/* error exit if file can't be opened                                */
/*********************************************************************/

FILE *open_file( char *name )
   {
    FILE *f;

    if ( strcmp( name, "-" ) == 0 )
        return( stdin );
    else
        if ( (f = fopen( name, "r" ) ) )
            return( f );
        else
           {
            perror( name );
            exit( 1 );
           }
   }

/* comp_char[c] contains the complementary nucleotide character for c */

char comp_char[] = {
                  /*          0     1     2     3     4     5     6     7   */
                  /* 000 */  'X',  'X',  'X',  'X',  'X',  'X',  'X',  'X',
                  /* 010 */  'X',  'X',  'X',  'X',  'X',  'X',  'X',  'X',
                  /* 020 */  'X',  'X',  'X',  'X',  'X',  'X',  'X',  'X',
                  /* 030 */  'X',  'X',  'X',  'X',  'X',  'X',  'X',  'X',
                  /* 040 */  'X',  'X',  'X',  'X',  'X',  'X',  'X',  'X',
                  /* 050 */  'X',  'X',  'X',  'X',  'X',  'X',  'X',  'X',
                  /* 060 */  'X',  'X',  'X',  'X',  'X',  'X',  'X',  'X',
                  /* 070 */  'X',  'X',  'X',  'X',  'X',  'X',  'X',  'X',
                  /* 100 */  'X',  'T',  'V',  'G',  'H',  'X',  'X',  'C',
                  /* 110 */  'D',  'X',  'X',  'M',  'X',  'K',  'N',  'X',
                  /* 120 */  'X',  'X',  'Y',  'S',  'A',  'X',  'B',  'W',
                  /* 130 */  'X',  'R',  'X',  'X',  'X',  'X',  'X',  'X',
                  /* 140 */  'X',  't',  'v',  'g',  'h',  'X',  'X',  'c',
                  /* 150 */  'd',  'X',  'X',  'm',  'X',  'k',  'n',  'X',
                  /* 160 */  'X',  'X',  'y',  's',  'a',  'X',  'b',  'w',
                  /* 170 */  'X',  'r',  'X',  'X',  'X',  'X',  'X',  'X',
                 };

/* close_file() closes the file unless its stdin */

void  close_file( FILE *f )
   {
    if ( f != stdin )
        fclose( f );
   }

char  ftag[MAXTAG] = "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN";
char  rtag[MAXTAG] = "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN";
int   num_ns;

void  init( void )
   {
    tag_len = tag1_len + tag2_len; 

    ftag[tag_len] = '\0';
    rtag[tag_len] = '\0';
    num_ns = tag_len;
   }

void  insert( int c )
   {
    char x = toupper( c );

    if ( ftag[0] == 'N' ) 
        num_ns--;
    if ( x == 'N' )
        num_ns++;

    memmove( ftag, ftag+1, tag_len-1 );   /* left shift forward tag buffer */
    ftag[tag_len - 1] = x;                /* and put the new guy on the end */

    memmove( rtag+1, rtag, tag_len-1 );
    rtag[0] = comp_char[x];    
   }


int   has_site( char *tag )
   {
    if ( strncmp( tag + mlyI_pos, MLYI_SEQ, 5 ) == 0 )
        return( 1 );
    else
        return( 0 );
   }

void  print_str( char *s, int n )
   {
    while ( n-- > 0 && *s != '\0' )
        putchar( *s++ );
   }

void  print_tag( char *tag, char dir )
   {
    print_str( tag, tag1_len );
    printf( " " );
    print_str( tag + tag1_len, tag2_len );
    printf( " %c %8d", dir, pos );
    if ( print_hdrs ) printf( " %s", hdr );
    putchar( '\n' );
   }


void  check_tag()
   {
    if ( num_ns == 0 )
       {
        if ( has_site( ftag ) )
            print_tag( ftag, 'f' );
        if ( has_site( rtag ) )
            print_tag( rtag, 'r' );
       }
   }

                               /****************/
                               /* Main Program */
                               /****************/

main( int argc, char **argv )
   {
    int          nfiles;
    char       **filenames;
    int          i;
    FILE        *f;
    int          c;
    char         last_char;
    char        *h;

    parse_args( argc, argv, &nfiles, &filenames );
    if ( verbose )
        printf( "hello tag len is %d, %d files\n", tag_len, nfiles );

    init();                             /* set up tag buffs */

    for ( i = 0; i < nfiles; i++ )      /* for each file */
       {
        if ( verbose )
            printf( "File: %s\n", filenames[i] );
        f = open_file( filenames[i] );                  /* open it */
        last_char = '\n';  
        while ( (c = fgetc( f )) != EOF )               /* read chars */
            if ( c == '>' && last_char == '\n' )        /* fasta header? */
               {                                        /* if so, reset */
                seq_num++;                              /* counters and skip*/
                pos = 0;
                if ( print_hdrs ) h = hdr;
                while ( (c = fgetc( f )) != EOF && c != '\n' )
                   {
                    if ( verbose ) putchar( c );
                    if ( print_hdrs ) *h++ = c;         /* size check! */
                   }
                if ( verbose ) putchar( '\n' );
                if ( print_hdrs ) *h = '\0';
                last_char = '\n';
               }
            else 
               {                                        /* otherwise if its */
                if ( isalpha( c ) )                     /* sequence, then */
                   {                                    /* insert into buffs*/
                    insert( c );                        /* and print tags */
                    check_tag();
                    pos++;
                   }
                last_char = c;
               }
        close_file( f );                  /* done this file */
       }
   }
