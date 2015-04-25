/* Program:     tagsearch.c                                                  */
/* Programmer:  Sean R. McCorkle                                             */
/*              Biology Dept., Brookhaven National Laboratory                */
/* Language:    C                                                            */
/*                                                                           */
/* Description: read tags (n-nmers) from a file and then scans sequences     */
/*              (FASTA format) and reports occurrances.                      */
/*                                                                           */
/* Usage:       tagsearch [-fhvV] [-p<str>] <tagfile> [<seq file> ...]       */
/*               where <tagfile> has one n-mer per line                      */
/*               if no sequence files are specified, stdin is scanned.       */
/*               - can also be used to indicate stdin                        */
/*                                                                           */
/* $Id: tagsearch.c,v 0.2 2003/12/12 15:48:09 mccorkle Exp mccorkle $        */
/*                                                                           */
/*****************************************************************************/

#include <ctype.h>
#include <errno.h>
#include <libgen.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

static char tag_rcs_id[] =
    "$Id: tagsearch.c,v 0.2 2003/12/12 15:48:09 mccorkle Exp mccorkle $";

#define  MAX_HDR_LEN  256
#define  MAX_TAG_LEN  256
#define  MAX_STR_LEN  256

#ifndef NULL
#define NULL        ((void *) 0)
#endif

char  tagfile[MAX_STR_LEN+1];
int   tag_len;

int    verbose = 0;          /* set by -v */
char  *prefix = "";          /* set by -p */
int    print_filenames = 0;  /* set by -f */

void  version( void )
   {
    char *v;
    v = strdup("$Revision: 0.2 $"+11 );
    v[strlen(v)-2] = '\0';
    printf( "tagsearch v%s\n", v );
   }

void  help( void )
   {
    version();
    printf( "no help\n");
   } 


void  parse_args( int argc, char **argv, char *tagfile,
                  int *nfiles, char ***filenames )
   {
    extern char *optarg;
    extern int   optind;
    int          c;
    static char *def_files[] = { "-", "" };

    while ( (c = getopt( argc, argv, "fhp:vV" ) ) != -1 )
        switch( c )
           {
            case  'f':   print_filenames = 1;
                         break;
            case  'h':   help();
                         exit( 0 );
            case  'p':   prefix = strdup( optarg );
                         break;
            case  'v':   verbose = 1;
                         break;
            case  'V':   version();
                         exit( 0 );
            default:     help();
                         exit( 1 );
           }
    argc -= optind;
    argv += optind;
    if ( argc <= 0 )
       {
        fprintf( stderr, "no tag file specified.  (prosearch -h for help)\n" );
        exit( 1 );
       }
    /* Note to self:  needed here is a length check and a check on allowable */
    /* characters in pattern                                                 */
    strncpy( tagfile, *argv, MAX_STR_LEN );
    argc--;
    argv++;
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

void  rc( char *r, char *s )
   {
    char *p;

    p = s + strlen( s );
    while ( p > s )
       *r++ = comp_char[*(--p)];
    *r = '\0';
   }

                          /********************/
                          /* Permutation tree */
                          /********************/

int ind[128];  /* maps a,g,c,t chars to 0,1,2,3 indices, set in init_tree() */

char nt[] = { 'A', 'C', 'G', 'T' };  /* reverse map of ind[] */

typedef struct tnode {
                       struct tnode *ch[4];
                       char         *desc;
                     } TNODE;

TNODE *t_root;            /* tree root */

TNODE *new_tnode( void )
   {
    TNODE *t;

    if ( !( t = (TNODE *) malloc( sizeof(TNODE) ) ) )
       {
        perror( "can't allocate a tree node" );
        exit( errno );
       }
    t->ch[0] = t->ch[1] = t->ch[2] = t->ch[3] = NULL;
    t->desc = NULL;
    return( t );
   }

/*******************/
/* initialize tree */
/*******************/

void  init_tree( void )
   {
    int i;

    t_root = new_tnode();

    for ( i = 0; i < 128; i++ )
        ind[i] = -1;   
    ind['a'] = ind['A'] = 0;
    ind['c'] = ind['C'] = 1;
    ind['g'] = ind['G'] = 2;
    ind['t'] = ind['T'] = 3;
    if ( verbose )
       printf( "tree created\n" );
   }


/************************************************************************/
/* search tree for string.  If found, return value string for the match */
/* match, otherwise return NULL                                         */
/************************************************************************/

char *tree_lookup( char *s )
   {
    TNODE *t;

    if ( verbose )
        printf( "lookup [%s]\n", s );
    t = t_root;
    while ( *s != '\0' && ind[*s] >= 0 && t->ch[ind[*s]] !=NULL )
        t = t->ch[ind[*s++]];

    if ( *s == '\0' )
        return( t->desc );
    else if ( ind[*s] < 0 )
        return( NULL );       /* 'X' or anything else in s fails */
    else
        return( NULL );
   }

/********************************************/
/* make one entry into the permutation tree */
/********************************************/

void  tree_enter( TNODE *t, char *s, char *desc )
   {
    if ( verbose ) printf( "insert [%s]\n", s );
    if ( *s == '\0' )
       {
        t->desc = strdup( desc );
        /* printf( "inserted: %s\n", desc );  */
       }
    else
       {
        if ( t->ch[ind[*s]] == NULL )
           {
            /* printf( "need new node\n" );  */
            t->ch[ind[*s]] = new_tnode();
           }
        tree_enter( t->ch[ind[*s]], s+1, desc );
       }
   }

/****************************************************************************/
/* put a new tag entry and its reverse-complement into the permutation tree */
/****************************************************************************/

void  enter_tag( char *ftag )
   {
    static  char  rtag[MAX_TAG_LEN+1];
    static  char  desc[MAX_STR_LEN+1];

    if ( verbose ) printf( "enter table [%s]\n", ftag );
    sprintf( desc, "f %s", ftag );
    tree_enter( t_root, ftag, desc );

    rc( rtag, ftag );
    sprintf( desc, "r %s", ftag );
    tree_enter( t_root, rtag, desc );
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
            exit( errno );
           }
   }

/* close_file() closes the file unless its stdin */

void  close_file( FILE *f )
   {
    if ( f != stdin )
        fclose( f );
   }

void  load_tags( char *tagfile )
   {
    FILE        *f;
    static char  line[MAX_STR_LEN+1];
    static char  tag[MAX_TAG_LEN+1];
    int          len;
    static char  fulltag[2*MAX_TAG_LEN+1];

    if ( verbose )
        printf( "tagfile: %s\n", tagfile );
    if ( verbose ) 
        printf( "tag_len is %d\n", tag_len );

    init_tree();

    f = open_file( tagfile );
    tag_len = -1;
    while ( fgets( line, MAX_STR_LEN, f ) )
       {
        if ( sscanf( line, "%s", tag ) > 0 )
           {
            len = strlen( tag );
            if ( verbose )
                printf( "tag [%s] %d\n", tag, len );
            if ( tag_len < 0 )
                tag_len = len;
            else if ( tag_len != len )
               {
                fprintf( stderr,
                         "file: %s, tag \"%s\" length %d differs from %d\n", 
                                tagfile, tag, len, tag_len );
                exit( 1 );
               }
            strcpy( fulltag, prefix );   /* note: need to put protections*/
            strcat( fulltag, tag );
            enter_tag( fulltag );
           }
       }
    close_file( f );
    tag_len += strlen( prefix );
   }

                       /**********************/
                       /* comparision buffer */
                       /**********************/

static char  buff[MAX_TAG_LEN+1];

void  init_buff( void )
   {
    size_t i;
    for ( i = 0; i < tag_len; i++ )
        buff[i] = 'X';
    buff[tag_len] = '\0';
   }

void  enter_base( int c )
   {
    size_t i;

    for ( i = 1; i < tag_len; i++ )   /* memmove would be better */
        buff[i-1] = buff[i];
    buff[tag_len-1] = c;
   }


/*****************************************************/
/* scan_file() - open file and process its sequences */
/*****************************************************/

void  scan_file( char *filename )
   {
    FILE        *f;
    int          c, last;
    static char  hdr[MAX_HDR_LEN+1];
    int          pos;   /* position within string */
    char        *rec;

    if ( verbose )
        printf( "file %s\n", filename );
    f = open_file( filename );
    last = '\n';
    while ( (c = fgetc( f ) ) != EOF )
       {
        if ( last == '\n' && c == '>' )
           {
            fgets( hdr, MAX_HDR_LEN, f );
            if ( verbose )
                printf( "seq: %s", hdr );
            init_buff();
            pos = -tag_len;
            last = '\n';
           }
        else if ( isalpha( c ) )
           {
            enter_base( c );
            ++pos;
            if ( rec = tree_lookup( buff ) ) 
               {
                printf( "%s %s %10d ",  buff, rec, pos );
                if ( print_filenames )
                    printf( " %s", basename( filename ) );
                printf( " %s", hdr );
               }
            last = c;
           }
        else
            last = c;
       }
    close_file( f );
   }



                              /****************/
                              /* Main Program */
                              /****************/


main( int argc, char **argv )
   {
    int     nfiles;
    char  **filenames;
    int     i;
    char    tagfile[MAX_STR_LEN+1];

    parse_args( argc, argv, tagfile, &nfiles, &filenames );

    load_tags( tagfile );

    for ( i = 0; i < nfiles; i++ )
        scan_file( filenames[i] );
   }
