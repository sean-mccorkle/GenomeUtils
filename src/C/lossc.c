/* Program:      lossc                                                       */
/* Programmers:  Lawrence Brenner, initial development Aug, 2001.            */
/*               updated by Sean McCorkle, Mar 2002                          */
/* Language:     C                                                           */
/*                                                                           */
/* Description:  Lots of short sequence comparisons.                         */
/*                                                                           */
/* Usage:        lossc [-hmV] [-v<n>] <thresh> <data seqs> [<db seqs>]       */
/*                                                                           */
/*               <thresh> is an integer indicated the maximum edit distance  */
/*                        considered for a match                             */
/*               <data seqs> is a file of short sequences and descriptions   */
/*               <db seqs>  is a file of short sequences and descriptions    */
/*                                                                           */
/*                 -f<c>  change desc seperator to <c> (default is '|')      */
/*                 -i<n>  use <n> for indel penalty (default is 1)           */
/*                 -m     show misses as well as hits                        */
/*                 -s<n>  use <n> for substitution penalty (default is 1)    */
/*                 -v<n>  verbosity level n (1, 2, ...) higher means more    */
/*                 -h     print usage message, then exit                     */
/*                 -V     print version, then exit                           */
/*                                                                           */
/* $Id: lossc.c,v 1.5 2003/01/19 13:35:11 mccorkle Exp mccorkle $ */
/*****************************************************************************/

#include <errno.h>
#include <libgen.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifndef SOLARIS        /* Sun's are really a pain in the neck sometimes */
#include <unistd.h>
#endif


static char const rcsid[] =
    "$Id: lossc.c,v 1.5 2003/01/19 13:35:11 mccorkle Exp mccorkle $";
static char rcsvers[] = "$Revision: 1.5 $";


#define   MAX_SEQ_SIZE       100     /* max sequence length (matrix size) */
#define   LARGE             (MAX_SEQ_SIZE + 1)
#define   MAX_FILENAME       255     /* max length of filename */
#define   MAX_MSG_LEN        255     /* max length of filename */
#define   MAX_NUM_SEQS    200000     /* maximum number of sequences */
#define   MAX_DESCRIP_LEN   1000     /* cutoff for descriptions */
#define   MAX_LINE_LEN      (MAX_SEQ_SIZE + MAX_DESCRIP_LEN + 10)

typedef struct {  char *seq;         /* pointer to tag sequence */
                  char *descrip;     /* pointer to tag description */
               } TAGREC;

int         min_seq_length;           /* set by main() */
char        field_sep = '|';          /* set by -f<char> */
int         misses;                   /* set by -m option in get_args */
int         verbosity;                /* set by -v option in get_args */
int         indel_penalty = 1;        /* set by -i option in get_args */
int         subst_penalty = 1;        /* set by -s option in get_args */
char       *progname = "lossc";       /* set by get_args */

static int  ed_dist[MAX_SEQ_SIZE][MAX_SEQ_SIZE]; /*The Smith-Waterman matrix*/


/**********************************************/
/* usage() - print usage message and bomb out */
/**********************************************/

void  usage( void )
   {
    printf( "Usage:  lossc [-hmV] [-f<c>] [-i<n>] [-s<n>] [-v<n>] <thresh> <seqfile> [<seqfile2>]\n" );
    exit( 1 );
   }

void  *malloc_safely( size_t size )
   {
    void        *p;

    if ( ! (p = malloc( size )) )
       {
        fprintf( stderr, "%s: failed to malloc %d bytes\n(run again with -v2 option for estimating shortfall)", progname, size );
        exit( 1 );
       }
    return( p );
   }

/*************************************************************************/
/* allocates storage (safely) for a new copy of string s, or a substring */
/* of s of length n.  If n < 0, the whole string is safely duplicated.   */
/*************************************************************************/

char  *dup_safely( char *s, int n )
   {
    char *t;

    if ( n < 0 )
        n = strlen( s );
    t = (char *) malloc_safely( n + 1 );
    strncpy( t, s, n );
    t[n] = '\0';            /*  need to terminate strncpy() */
    return( t );  
   }

/**************************************************/
/* repeat_char( c, n ) - return a string of n c's */
/**************************************************/

char *repeat_char( char c, int n )
   {
    char *s;
    int   i;

    s = malloc_safely( n + 1 );
    for ( i = 0; i < n; i++ )
        s[i] = c;
    s[i] = '\0';
    return( s );
   }

/**************************************************************************/
/* get_args() - parse args for options and set global values accordingly, */
/* and then extract edit_threshold and filenames and return them.         */
/**************************************************************************/

void get_args( int     argc,         /* input- argc from main( argc, argv )*/
               char   *argv[],       /* input- argv from main( argc, argv )*/
               int    *edit_thresh,  /* output- edit distance threshold    */
               char   *file1,       /* output- mandatory filename         */
               char   *file2        /* output- opt. 2nd file or "" if none*/
            )
   {
    extern char *optarg;
    extern int   optind;
    int          c;

    progname = strdup( basename( *argv ) );
    misses = 0;
    verbosity = 0;
    while ( (c = getopt( argc, argv, "f:hi:ms:v:V")) != -1 )
        switch ( c )
           {
            case 'f':  field_sep = *optarg;
                       break;
            case 'i':  indel_penalty = atoi( optarg );
                       break;
            case 'm':  misses = 1;
                       break;
            case 's':  subst_penalty = atoi( optarg );
                       break;
            case 'v':  verbosity = atoi( optarg );
                       break;
            case 'V':  
                       rcsvers[strlen(rcsvers)-1] = '\0';
                       printf( "lossc, %s\n", rcsvers+1 );
                       exit( 0 );
            case 'h':  
            default:
                       usage();
           }
    argc -= optind;
    argv += optind;

    if ( argc < 2 )
        usage();

    *edit_thresh = atoi( argv[0] );
    strncpy( file1, argv[1], MAX_FILENAME );  /* note: add termination here */
    if ( argc > 2 )
        strncpy( file2, argv[2], MAX_FILENAME );
    else
        *file2 = '\0';

    if ( verbosity > 1 )
        printf( "file1 [%s] file2 [%s]\n", file1, file2 );
   }


/* note to self: put in lots of input checks, including max seq length */
/* check. */

void  load_sequences(  char   *filename, /* input- short seq filename        */
                       TAGREC *tags,     /* output- array of seqs & descrips */
                       int    *count     /* output- number of tags in file   */
                    )
 {
  FILE          *in;                     /* input file handle */
  static char    line[MAX_LINE_LEN+1];   /* input line buffer */
  int            len;                    /* input line length */
  int            line_num = 0;           /* line number for diagnostic output*/
  char          *d;
  int            k;
  static char    msg[MAX_MSG_LEN+1];

  sprintf( msg, "(%s) %s ", progname, filename );  /* check on sizes */
  if ( ! (in = fopen( filename, "r" ) ) )
    {
     perror( msg );
     exit( errno );
    }

  *count = 0;
  while ( ( *count < MAX_NUM_SEQS ) && fgets( line, MAX_LINE_LEN, in ) )
    {
     line_num++;
     len = strlen( line );
     if ( ( len > 0 ) && ( line[len-1] == '\n' ) ) /* chomp ending '\n' */
         line[--len] = '\0';
      d = strchr( line, ' ' );  /* note: update this to include tabs too */
      k = (int) ( d - line );
      if ( verbosity > 1 )
          printf( "k %d d %x [%s] %d\n", k, d, line, len );
      if ( k > MAX_SEQ_SIZE )
         {
          fprintf( stderr, "(%s) %s, line %d: sequence length > %d\n%s\n",
                           progname, filename, line_num, MAX_SEQ_SIZE, line );
          exit( 1 );
         }
      if ( k > 0 )
        {
         tags->seq = dup_safely( line, k );
         tags->descrip = dup_safely( d + 1, -1 );
         if ( verbosity > 1 )
   	     printf( "input: [%s] [%s]\n", tags->seq, tags->descrip );
         tags++, (*count)++;
	}
      else if ( len > 0 )
        {
         tags->seq = dup_safely( line, -1 );
         tags->descrip = " ";
         if ( verbosity > 1 )
   	     printf( "input: [%s] [%s]\n", tags->seq, tags->descrip );
         tags++, (*count)++;
        }
    }
    fclose( in );
    if ( *count >= MAX_NUM_SEQS )
       {
        fprintf( stderr, "(%s) file %s exceeds limit of %d sequences\n",
                           progname, filename, MAX_NUM_SEQS );
        exit( 1 );
       }
    if ( verbosity > 0 )
        printf( "%d sequences in file %s\n", *count, filename );
} 

/*******************************************************************/
/* min_length() returns the minium length of the seqs in the array */
/* note that this expects the initial value of min as an argument  */
/*******************************************************************/

int   min_length( TAGREC *tags, int n, int min )
   {
    int  len;

    while ( n-- > 0 )
       {
        len = strlen( (tags++)->seq );
        if ( len < min )
            min = len;
       }
    if ( verbosity > 0 )
        printf(  "Min seq length %d\n", min );
    return( min );
   }

/*******************************************************/
/* trunc_lengths() truncates all seqs down to size len */
/*******************************************************/

void  trunc_lengths( TAGREC *tags, int n, int len )
   {
    char *seq;

    while ( n-- > 0 )
       {
        seq = (tags++)->seq;
        if ( strlen( seq ) > len )
           {
            if ( verbosity > 1 ) printf( "    truncating [%s] ", seq );
            seq[len] = '\0';
	    if ( verbosity > 1 ) printf( " to [%s]\n", seq );
           }
       }
   }

               /*****************************/
               /* Edit distance calculation */
               /*****************************/

/*****************************************************************/
/* Fills the very top row and left column with increasing scores */
/*****************************************************************/

void  initialize_matrix( void )
   {
    int row;
    int col;

    for ( row = 0; row < MAX_SEQ_SIZE; row++ ) /*Sets the values of r,0 */
        ed_dist[row][0] = row;

    for ( col = 1; col < MAX_SEQ_SIZE; col++ ) /*Sets the values of 0,c */
        ed_dist[0][col] = col;
   }


/***************************************************************************/
/* update_min() - selects the minimum score at ed_dist[row,col], which may */
/* necessitate checking seq1[row-1] and seq2[col-1].  Returns the minimum  */
/* value in final                                                          */
/***************************************************************************/

void  update_min( int     row,   /* input- row number for evaluation */
                  int     col,   /* input- col number for evaluation */
                  char  *seq1,   /* input- 1st short sequence        */
                  char  *seq2,   /* input- 2nd short sequence        */
                  int   *final    /* output- minimum value            */
                )
   {
    int           b = 0;          /*Used in the substitution, match function.*/
  

    if ( seq1[row - 1] == seq2[col - 1] ) /*Checks for matches. */
        b = 0;
    else                                     /* Checks for substitutions. */
        b = subst_penalty;
  
    *final = ed_dist[row - 1][col - 1] + b;    /* - 1; */
  
    if ( ed_dist[row - 1][col] + 1 < *final )  /* Checks for insertions. */
        *final = ed_dist[row - 1][col] + indel_penalty;
    if ( ed_dist[row][col - 1]  + 1 < *final ) /* Checks for deletions. */
        *final  = ed_dist[row][col - 1] + indel_penalty; 

    ed_dist[row][col] = *final; /*++final; */
   }

/*****************************************************************************/
/* This version of edit_distance counts on the previous values of            */
/* ed_dist[r][c] from the previous sequence comparsion to be present,        */
/* up to the subsquare determined by start_diag.                             */
/*                                                                           */
/* returns:  dist - edit distance (or threshold + 1 if aborted)              */
/*           last_diag - matrix diag indicates where calculation stopped     */
/*                       (< seq.length() if aborted early because threshold  */
/*                        exceeded)                                          */
/*****************************************************************************/

void edit_distance( char *seq1,       /* input- first short sequence         */
		    char *seq2,       /* input- 2nd short sequence           */
		    int   threshold,  /* input- abort if dist > threshold    */
		    int   start_diag, /* input- recalculate mat @ this column*/
		    int  *last_diag,  /* output- last diag calc. before abort*/
		    int  *dist )      /* output- resultant distance          */
   {
    int  eta = 0;      /* Minimum edit distance changes.                     */
    int  min;          /* The lowest calculated value at one position        */
    int  r, c;         /* row and column counters                            */
 
    int  rows = min_seq_length;  /* number of rows in the matrix             */
    int  cols = min_seq_length;  /* number of columns                        */
    int  diag;         /* Variable that analyzes only the revelant diagonals.*/

    /*************************************************************************/
    /* Begin filling the matrix ed_dist[r][c] in subsquares, starting        */
    /* in the upper right corner, or actually at the specified start position*/
    /*************************************************************************/

    for ( diag = start_diag; ((diag <= rows) && (eta <= threshold)); diag++ )
       {
        eta = ed_dist[0][diag];
        for (  r = 1;  r < diag;  r++  )         /* do the column */
	   {
	    update_min( r, diag, seq1, seq2, &min );
	    if ( min < eta )
	        eta = min;
           }

        for (  c = 1;  c <= diag;  c++  )        /* now do the row */
	   {
	    update_min( diag, c, seq1, seq2, &min );
            if ( min < eta )
               eta = min;
           }
       }
 
    if ( eta > threshold )
        *dist = threshold + 1;    
    else
        *dist = ed_dist[rows][cols];

    *last_diag = diag;
   }


/****************************************************************************/
/* first_diff( a, b ) returns the index of the first position where         */
/* strings a & b differ.  (Note: we assume both strings are the same length */
/* - is this always going to be true?                                       */
/****************************************************************************/

int  first_diff( char *a, char *b )
   {
    int i = 0; 

    while ( ( i < min_seq_length ) && ( a[i] == b[i] ) )  
        i++;
    return i;
   }


void  cross_compare( int      edit_thresh,
                     TAGREC  *tags, 
                     int      n_tags, 
                     TAGREC  *db,
                     int      n_db
                   )
   {
    int             t;             /* database sequence index                */
    int             d;             /* counter for the database.              */
    int             ed;            /* edit distance between the sequences.   */
    int             same;          /* difference between the two sequences.  */
    int             last = LARGE;  /* last col position calculated in matrix.*/
    int             best_ed;
    char           *spacer;

    spacer = repeat_char( '-', min_seq_length );
    for ( t = 0; t < n_tags; t++ ) /* Checks the tags to the database. */
       {
        best_ed = edit_thresh + 1;
        for ( d = 0; d < n_db; d++ )
           { 
            if ( d == 0 )
                same = 0;
            else
                same = first_diff( db[d - 1].seq, db[d].seq ); 

            /***************************************************************/
            /* if we aborted early the previous time, and this sequence    */
            /* is the same BEYOND that point, then there's no need to      */
            /* to bother - this one doesn't match in the threshold either. */
            /***************************************************************/

            if ( same < last )
               {
                edit_distance( tags[t].seq, db[d].seq, edit_thresh, same, 
                               &last, &ed );
                if ( ed <= edit_thresh )
                   {
                    printf( "%s %s %d %s%c%s\n", tags[t].seq, db[d].seq, ed,
                               tags[t].descrip, field_sep, db[d].descrip );
                    best_ed = ed;
                   }
               }
           }
        if ( misses && ( best_ed > edit_thresh ) )  /* no hits, print miss */
            printf( "%s %s X %s\n", tags[t].seq, spacer, tags[t].descrip );

       }   
   }


void  self_compare( int      edit_thresh,
                    TAGREC  *tags,
                    int      n_tags )
   {
    int             i, j;
    int             same;          /* difference between the two sequences.  */
    int             last = LARGE;  /* last col position calculated in matrix.*/
    int             ed;
    int             best_ed;
    char           *spacer;

    spacer = repeat_char( '-', min_seq_length );
    if ( verbosity > 1 )
        printf( "self compare\n" );

    for ( i = 0; i < n_tags; i++ )
       {
        best_ed = edit_thresh + 1;
        for ( j = i+1; j < n_tags; j++ )
           {
            if ( j == i + 1 )
                same = 0;
            else
                same = first_diff( tags[j - 1].seq, tags[j].seq ); 

            /***************************************************************/
            /* if we aborted early the previous time, and this sequence    */
            /* is the same BEYOND that point, then there's no need to      */
            /* to bother - this one doesn't match in the threshold either. */
            /***************************************************************/

            if ( same < last )
               {
                edit_distance( tags[i].seq, tags[j].seq, edit_thresh, same, 
                               &last, &ed );
                if ( ed <= edit_thresh )
                   {
                    printf( "%s %s %d %s%c%s\n", tags[i].seq, tags[j].seq, ed,
                               tags[i].descrip, field_sep, tags[j].descrip );
                    best_ed = ed;
                   }
               }
           }
        if ( misses && ( best_ed > edit_thresh ) )  /* no hits, print miss */
            printf( "%s %s X %s\n", tags[i].seq, spacer, tags[i].descrip );
       }
   }


void  dump_tags( TAGREC *tags, int n )
   {
    int  i;
    for ( i = 0; i < n; i++ )
        printf( "%8d %s %s\n", i, tags[i].seq, tags[i].descrip );
   }

int  tagseq_cmp( const void *a, const void *b )
   { 
    return( strcmp( ((TAGREC *)a)->seq, ((TAGREC *)b)->seq ) );
   }

void  sort( TAGREC *tags, int n )
   {
    if ( verbosity > 0 )
       {
        printf( "sorting %d tags...\n", n );
        if ( verbosity > 1 )  dump_tags ( tags, n );
       }

    qsort( (void *) tags, n, sizeof( TAGREC ), tagseq_cmp );

    if ( verbosity > 0 )
       {
        printf( "...done\n" );
        if ( verbosity > 1 )  dump_tags ( tags, n );
       }
   }

                           /****************/
                           /* Main Program */
                           /****************/


int  main( int argc, char *argv[] )
  
   {
    static char   tag_file[MAX_FILENAME+1]; /* tag file name                 */
    static TAGREC tags[MAX_NUM_SEQS];       /* array of tags to be checked   */
    int           n_tags;                   /* number of tags in array       */

    static char   db_file[MAX_FILENAME+1];  /* optional tag database filename*/
    static TAGREC db[MAX_NUM_SEQS];         /* array of database tags        */
    int           n_db;                     /* number tags in database array */

    int           edit_thresh;     /* How much of an edit distance you want. */

    get_args( argc, argv, &edit_thresh, tag_file, db_file );

    initialize_matrix();

    load_sequences( tag_file, tags, &n_tags ); /* read data tags*/
    min_seq_length = min_length( tags, n_tags, LARGE );
    if ( *db_file != '\0' )                    /* read optional database tags*/
       {
        load_sequences( db_file, db, &n_db ); 
        min_seq_length = min_length( db, n_db, min_seq_length );
        trunc_lengths( tags, n_tags, min_seq_length );        
        trunc_lengths( db, n_db, min_seq_length );        
        sort( db, n_db );                      /* sort alphabetically */
        cross_compare( edit_thresh, tags, n_tags, db, n_db );
       }
    else
       {
        trunc_lengths( tags, n_tags, min_seq_length );        
        sort( tags, n_tags );
        self_compare( edit_thresh, tags, n_tags );
       }

   }



