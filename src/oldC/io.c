/* Module:       io.c                                                */
/* Programmer:   Sean R. McCorkle                                    */
/* Language:     C                                                   */
/* Description:  Routines for handling sequence I/O.                 */
/*                                                                   */
/* $Id: io.c,v 3.9 2000/03/16 19:48:30 mccorkle Exp mccorkle $       */
/*********************************************************************/


static char io_rcs_id[] =
      "$Id: io.c,v 3.9 2000/03/16 19:48:30 mccorkle Exp mccorkle $";

#include <ctype.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "seqlib.h"

#define  MAX_MSG_LENGTH  200

void  make_char_set_index( char *str, char *char_set );
void  get_seq_buff( int s );
int   get_next_char_from_seqbuff( int s );
void  tab( FILE *output, int n );

/* array (hopefully small!) of buffers, flags, etc for handling multiple */
/* open sequence files.                                                  */

#define MAX_OPEN_SEQF  5

char   *seq_buff[MAX_OPEN_SEQF];
int     at_eof[MAX_OPEN_SEQF];
char   *next_p[MAX_OPEN_SEQF];
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
        errno = ENOSR;
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
    seq_buff[s] = (char *) malloc_safely( MAX_LINE_LENGTH, 
                                          "occurred in open_seqf, seq_buff" );
    *seq_buff[s] = '\0';              /* indicates "no data yet"      */
    next_p[s] = seq_buff[s];          /* also indicates "no data yet  */
    make_char_set_index( allowed, is_allowed[s] );
    make_char_set_index( ignored, is_ignored[s] );
    in_use[s] = 1;
    at_eof[s] = 0;
    lineno[s] = 0;
    return( s );
   }

/****************************************************/
/* close a sequence channel (freeing it for re-use) */
/****************************************************/

int  close_seqf( int s )
   {
    if ( s >= 0 && s < MAX_OPEN_SEQF )
       {
        free( seq_buff[s] );
        free( seqfilename[s] );
        in_use[s] = 0;
        return( 1 );
       }
    else
        return( 0 );
   }

/* call free_lines() to free all the storage allocated for the lines */
/* *nlines = 0 if no header was read                                 */

void load_header( int s, int *nlines, char ***lines )

   {
    static   char *loclines[MAX_HEADER_LINES];
    int      n = 0;

    if ( ! *seq_buff[s] )          /* if theres nothing in buffer, read */
        get_seq_buff( s );         /* first line                        */

    while ( (! at_eof[s] ) && n < MAX_HEADER_LINES && 
         ( *seq_buff[s] == '>' || *seq_buff[s] == ';' || blank( seq_buff[s] )))
       {
        if ( *seq_buff[s] == '>' || *seq_buff[s] == ';' )
            loclines[n++] = strdup( seq_buff[s] );
        get_seq_buff( s );                      /* get new line */
       }
    if ( n >= MAX_HEADER_LINES )
       {
        fprintf( stderr, 
                 "%s (load_header): file %s header exceeds %d lines\n", 
                 progname, seqfilename[s], MAX_HEADER_LINES );
        exit( 1 );
       }
    *nlines = n;
    if ( n < 1 )
        *lines = (char **) 0;
    else
       {
        *lines = (char **) malloc_safely( n * sizeof( char * ),
                                          "occurred in load_header" );
        while ( n-- > 0 )
            (*lines)[n] = loclines[n];
       }
   }


void skip_header( int s )

   {
    if ( ! *seq_buff[s] )       /* if theres nothing in buffer, read */
        get_seq_buff( s );      /* first line                        */

    while ( (! at_eof[s]) && 
       ( *seq_buff[s] == '>' || *seq_buff[s] == ';' || blank( seq_buff[s] ) ) )
        get_seq_buff( s );
   }

/* Frees all memory allocated for a lines structure, like that returned  */
/* by load_header()                                                      */

void free_lines( int nlines, char **lines )

   {
    while ( nlines-- > 0 )
        free( lines[nlines] );
    if ( lines )
        free( lines );
   }


/***************************************************************************/
/* load sequence - reads one of our modified FASTA-style ascii sequence    */
/* files, globbing in the entire sequence in one shot.                     */
/* Returns -1 if an error occured (sequence length exceeds max_seq_size)   */
/*          0 if no sequence was read (EOF) encountered                    */
/*          1 if a sequence was read successfully                          */
/* returns number of characters in sequence read (not counting ignored     */
/* characters.  0 means no characters read, -1 if sequence length exceeds  */
/* max_seq_size or if a character was read which is not in the allowed     */
/* set string, or the ignored set string                                   */
/***************************************************************************/

int load_seq( int   s,        /* input; input seq file desc (from open_seqf) */
              int   max_seq_size, /* input; size of seq array                */
              char *seq,      /* output; where to place output sequence      */
	      int  *seq_len   /* output; actual length read                  */
            )  

   {
    int        c;
    int        n;

    n = 0;
    while ( ( n < max_seq_size ) && 
            ( c = get_next_sequence_char( s ) ) &&
            c != ERROR && c != EOF )
       {
        n++;
        *seq++ = (char) c;
       }

    if ( n >= max_seq_size )
       {
        fprintf( stderr, "%s: input sequence in file \"%s\" too long! (>%d)\n",
                 progname, seqfilename[s], max_seq_size );
        return( -1 );
       }
    if ( c == ERROR )
        return( -1 );
    if ( (c == EOF || c == 0) && n == 0 )
        return( 0 );
    *seq = '\0';
    *seq_len = n;
    return( 1 );
   }

/* returns EOF, ERROR, '\0' or the next allowed character (removes ignoreds) */

int  get_next_sequence_char( int s )
   {
    int c;

    do {
        c = get_next_char_from_seqbuff( s );
       }
    while ( c && c != EOF && is_ignored[s][(char) c] );

    if ( ( ! c ) || c == EOF )
        return( c );

    if ( is_allowed[s][(char) c] )
        return( c );
    else
       {
        fprintf( stderr,
           "%s: bad character \"%c\" (%x hex) in sequence file %s, line %d\n", 
           progname, (char) c, c, seqfilename[s], lineno[s] );
        return( ERROR );
       }
   }

/* return the next raw character sitting in the buffer for channel s.  Reload*/
/* the buffer if necessary. returns EOF if end of file is hit on reload,     */
/* 0 if the buffer reload gets a next header ("> ") line, or the             */
/* actual raw character.  This skips past comment lines ("; ")               */

int  get_next_char_from_seqbuff( int s )
   {
    if ( at_eof[s] )
        return( EOF );

    if ( (! *next_p[s]) || next_p[s] > seq_buff[s] + MAX_LINE_LENGTH )
       {
        do {
            get_seq_buff( s );
           }
        while ( (! at_eof[s]) && *seq_buff[s] == ';' );

        if ( at_eof[s] )
            return( EOF );
        if ( *seq_buff[0] == '>' )
            return( '\0' );
       }
    return( *(next_p[s])++ );
   }


/*****************************************************************************/
/* get_seq_buff( s ) - read another line from file desc s (from open_seqf),  */
/* into seq_buff for file s.                                                 */
/* set global boolean at_eof[s] when end-of-file is encountered.             */
/* Note: this doesn't handle lines longer than MAX_LINE_LENGTH properly!     */
/*       FIX!  Put in error handling here!                                   */
/*****************************************************************************/

void  get_seq_buff( int s )

   {
    if ( at_eof[s] )
        return;

    if ( ! fgets( seq_buff[s], MAX_LINE_LENGTH, seqf[s] ) )
       {
#ifdef DEBUG
        printf( "get_seq_buff: eof, line_no %d\n", lineno[s] );
#endif
        at_eof[s] = 1;
        return;
       }
    next_p[s] = seq_buff[s];
    lineno[s]++;
#ifdef DEBUG
    printf( "get_seq_buff: line %d [%s] (%d)\n", 
              lineno[s], seq_buff[s], at_eof[s] );
#endif
   }


void print_lines( FILE *output, int nlines, char **lines )

   {
    int  i;
    for ( i = 0; i < nlines; i++ )
        fputs( lines[i], output );
   }


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
/* */

void  print_seq( FILE *output, char *seq, int start, int stop, int indent,
                 int linelen, int spaces, int xtralines
               )

   {
    int  n, i, j;

    n = strlen( seq );
    if ( start < 1 )
        start = 1;
    if ( stop < 1 || stop > n )
        stop = n;
    if ( linelen < 1 )
        linelen = 60;
#ifdef DEBUG
    printf( "start is %d\n", start );
#endif
    tab( output, indent );
    for ( j = 1, i = start - 1; i < stop; i++, j++ )
       {
        putc( seq[i], output );
        if ( (j % linelen) == 0 )
	   {
            putc( '\n', output );
            tab( output, indent );
	   }
        else if ( (j % 10) == 0 )
            putc( ' ', output );
       }
    if ( (j-1) % 60 )
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

/* Safely routines begin here. */

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
