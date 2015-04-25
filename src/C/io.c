/* Module:       io.c                                                */
/* Programmer:   Sean R. McCorkle                                    */
/* Language:     C                                                   */
/* Description:  Routines for handling sequence I/O.                 */
/*                                                                   */
/* $Id: io.c,v 4.3 2003/01/14 21:54:08 mccorkle Exp mccorkle $       */
/*********************************************************************/


static char io_rcs_id[] =
      "$Id: io.c,v 4.3 2003/01/14 21:54:08 mccorkle Exp mccorkle $";

#include <ctype.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "seqlib.h"

#define  MAX_MSG_LENGTH  200

void  make_char_set_index( char *str, char *char_set );
void  tab( FILE *output, int n );

/* array (hopefully small!) of buffers, flags, etc for handling multiple */
/* open sequence files.                                                  */

#define MAX_OPEN_SEQF  5

int     in_use[MAX_OPEN_SEQF] = { 0, 0, 0, 0, 0 };
FILE   *seqf[MAX_OPEN_SEQF];
char   *seqfilename[MAX_OPEN_SEQF];
char    lineno[MAX_OPEN_SEQF];
char    is_allowed[MAX_OPEN_SEQF][NUM_CHARS];
char    is_ignored[MAX_OPEN_SEQF][NUM_CHARS];

/*************************************************************************/
/* Open a sequence file for reading by load_sequence() and load_header.  */
/* At the same time, establish which chars are allowed and which ignored */
/* returns a integer file specifier which is passed as an argument to    */
/* those routines, or -1 if an error occurs.  If filename is "-", then   */
/* stdin is used.                                                        */
/* Let the user beware:  if allowed and ignored share a character, its   */
/* trouble!                                                              */
/*************************************************************************/

int  open_seqf( char *filename,
                char *allowed,  /* input; legal chars for seq, ie. "acgt" */
                char *ignored   /* input; chars to ignore, ie. "\n\t"     */
              )

   {
    int s = 0; 

    while ( s < MAX_OPEN_SEQF && in_use[s] )
        s++;
    if ( s >= MAX_OPEN_SEQF )
       {
        errno = EMFILE;
        return( -1 );
       }
    if ( strcmp( filename, "-" ) == 0 )
       {
        seqf[s] = stdin;
        seqfilename[s] = strdup( "- (stdin)" );
       }
    else
       {
        if ( ! ( seqf[s] = fopen( filename, "r" ) ) )
            return( -1 );
        seqfilename[s] = strdup( filename );
       }
    make_char_set_index( allowed, is_allowed[s] );
    make_char_set_index( ignored, is_ignored[s] );
    in_use[s] = 1;
    lineno[s] = 0;
    return( s );
   }

/***************************************************************************/
/* read_seq() reads one bnl-fasta format sequence from the file s          */
/* returns number of characters in sequence read (not counting ignored     */
/* characters.  0 means no characters read, -1 means EOF encountered */
/* -2 if sequence length exceeds  */
/* max_seq_size or if a character was read which is not in the allowed     */
/* set string, or the ignored set string                                   */
/***************************************************************************/

/*******************************************/
/* Note: Not yet handling ; comments!!!!!! */
/*******************************************/
/* And put in line number counting for error messages */

int  read_seq( int   s,             /* input; seq file desc (from open_seqf)*/
               int   max_hdr_size,  /* input; size of header buffer */
               char *hdr,           /* output; header buffer */
               int   max_seq_size,  /* input: size of sequence buffer */
               char *seq            /* output; sequence buffer */
             )
   {
    int   c, last_c;
    FILE *f;
    int   n;

    if ( ! in_use[s] )
       {
        fprintf( stderr, "read_seq() called on unused seq file #%d\n", s );
        return( -2 );
       }
    f = seqf[s];
    if ( feof( f ) )
       return( -1 );
    c = getc( f );
    if ( c != '>' )
       {
        fprintf( stderr, "read_seq() did not find header #%d\n", s );
        return( -2 );
       }
	                                      /* read the header first */
    n = max_hdr_size;
    while ( n-- > 0 && (c = getc( f )) != EOF && c != '\n' )
        *hdr++ = c;
    *hdr = '\0';
    if ( n <= 0 )
       {
        fprintf( stderr, 
                 "read_seq(): header length exceeds maximum of %d\n",
                 max_hdr_size );
        return( -2 );
       }
                                             /* now get the sequence */

    n = 0;
    last_c = '\n';
    while ( n < max_seq_size  && (c = getc( f )) != EOF && 
            ! ( last_c == '\n' && c == '>' )
          )
       {
        if ( ! is_ignored[s][(char) c] )
           {
            if ( is_allowed[s][(char) c] )
                *seq++ = c;
            else
               {
                fprintf( stderr,
                         "read_seq(): bad character \"%c\" (%x hex) in sequence file %s\n", 
                         (char) c, c, seqfilename[s] );
                return( -2 );
               }
            n++;
           }
        last_c = c;            
       }
    if ( n >= max_seq_size )
       {
        fprintf( stderr, 
                 "read_seq(): sequence length exceeds maximum of %d\n",
                 max_seq_size );
        return( -2 );
       }
    if ( c == '>' )
        ungetc( c, f ); 
    *seq = '\0';               /* don't forget to end the sequence string! */
    return( n );       
   }

/****************************************************/
/* close a sequence channel (freeing it for re-use) */
/****************************************************/

int  close_seqf( int s )
   {
    if ( s >= 0 && s < MAX_OPEN_SEQF )
       {
        free( seqfilename[s] );
        in_use[s] = 0;
        return( 1 );
       }
    else
        return( 0 );
   }


/****************************************************************************/
/* print_seq( f, seq, start, stop, indent, linelen, spaces, xtralines       */
/*                                                                          */
/*     f     - (input) file output goes to                                  */
/*                                                                          */
/*  Input Sequence Parameters                                               */
/*                                                                          */
/*  seq   - (input) sequence string                                         */
/*  start - (input) index of string to start printing at (default = 0)      */
/*  stop  - (input) index of last character in string to print (defaults    */
/*                  to full string)                                         */
/*                                                                          */
/*  Output Sequence Parameters                                              */
/*                                                                          */
/*  indent    - (input) spaces to indent at the beginning of each line      */
/*  linelen   - (input) characters per line of output (not including indent */
/*                      or spaces (0 defaults to 60)                        */
/*  spaces    - (input) num characters between spaces                       */
/*                      (0 or >linelen means no spaces at all)              */
/*  xtralines - (input) number of extra blank line                          */
/*                                                                          */
/****************************************************************************/

void  print_seq( FILE *output, char *seq, int start, int stop, int indent,
                 int linelen, int spaces, int xtralines
               )

   {
    int  n, i, j;

    n = strlen( seq );
    if ( start < 0 )
        start = 0;
    if ( stop < 1 || stop > n )
        stop = n;
    if ( linelen < 1 )
        linelen = 60;
#ifdef DEBUG
    printf( "start is %d\n", start );
#endif
    tab( output, indent );
    for ( j = 1, i = start; i < stop; i++, j++ )
       {
        putc( seq[i], output );
        if ( (j % linelen) == 0 )
	   {
            putc( '\n', output );
            tab( output, indent );
	   }
        else if ( spaces > 0 && (j % spaces) == 0 )
            putc( ' ', output );
       }
    if ( (j-1) % linelen )
        putc( '\n', output );
   }

void tab ( FILE *output, int n )
   {
    while ( n-- > 0 )
        putc( ' ', output );
   }

int count_char_occurrences( char *set, char *seq, int start, int stop )

   {
    static char  is_in_set[NUM_CHARS];
    char        *s, *s_end;
    int          n = 0, len;

    len = strlen( seq );
    if ( start < 1 )
        start = 1;
    s = seq + start - 1;
    if ( stop < 1 || stop > len )
        stop = len;
    s_end = seq + stop;
    make_char_set_index( set, is_in_set );    
    while ( *s && s < s_end )
        if ( is_in_set[ *s++ ] )
            n++;
    return( n );
   }

/* convert a string set into a character lookup table for fast  */
/* evaluation of whether or not a character is present in a set */
/* (If str is 0, set is empty) */

void make_char_set_index( char *str, char *char_set )

   {
    memset( char_set, 0, NUM_CHARS );
    if ( str )
        while ( *str )
            char_set[ *str++ ] = 1;
   }

                      /*******************************/
                      /* Safely routines begin here. */
                      /*******************************/

FILE *fopen_safely( char *name, char *mode )

   {
    FILE        *f;
    static char  msg[MAX_MSG_LENGTH];

    sprintf( msg, "(%s) %s", progname, name );
    if ( strcmp( "-", name ) == 0 )
        return(  (mode[0] == 'r') ? stdin : stdout );

    if ( !( f = fopen( name, mode ) ) )
       {
        perror( msg );
        exit( 1 );
       }
    return( f );
   }

void *calloc_safely( size_t nmemb, size_t size, char *errmsg )
   {
    void        *p;
    static char  msg[MAX_MSG_LENGTH];

    sprintf( msg, "(%s) - failed to calloc %d elements of %d bytes", 
                  progname, nmemb, size );
    if ( ! (p = calloc( nmemb, size )) )
       {
        perror( msg );
        if ( errmsg && *errmsg )
            fprintf( stderr, "%s\n",errmsg );
        exit( 1 );
       }
    return( p );
   }

void *malloc_safely( size_t size, char *errmsg )
   {
    void        *p;
    static char  msg[MAX_MSG_LENGTH];

    sprintf( msg, "(%s) - failed to malloc %d bytes", progname, size );
    if ( ! (p = malloc( size )) )
       {
        perror( msg );
        if ( errmsg && *errmsg )
            fprintf( stderr, "%s\n",errmsg );
        exit( 1 );
       }
    return( p );
   }
