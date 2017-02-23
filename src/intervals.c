/* Program:     intervals                                                    */
/* Programmer:  Sean R. McCorkle                                             */
/*              Biology Dept. Brookhaven National Laboratory                 */
/* Language:    C                                                            */
/*                                                                           */
/* Description: This a C version of the perl script, which itself was a      */
/*              generalized version of truncseq                              */
/*                                                                           */
/* Usage:       intervals -i n-m,[n-m,...]   [sequence-file]                 */
/*              intervals -f file  [sequence-file]                           */
/*                                                                           */
/*              There can be only one input sequence, and it must be in      */
/*              fasto format.  If no file is specified, or is "-", stdin is  */
/*              read.                                                        */
/*                                                                           */
/*              One of the two mutualy exclusive options -i or -q must be    */
/*              used to specify the desired intervals.                       */
/*                                                                           */
/*              The first sequence position is 0.  Intervals n-m start and   */
/*              include n and stop at m, but do not include m.  They may     */
/*              overlap.                                                     */
/*                                                                           */
/* Options:        -i n-m,...  comma-seperated intervals are specified on the*/
/*                          command line.  Each inteval is a n-m pair of     */
/*                          hyphen-separated integers                        */
/*                                                                           */
/*                 -f file    read intervals from file, each pair on one line*/
/*                            ('-' or space separated).                      */
/*                                                                           */
/* Compiling:   cc -O -o intervals intervals.c should do the job.            */
/*                                                                           */
/*****************************************************************************/


#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <errno.h>
#include <unistd.h>

static char intervals_rcs_id[] =
    "$Id: intervals.c,v 0.4 2007/11/26 16:46:56 mccorkle Exp mccorkle $"; 

#define MAX_NUM_INTS    1000000   /* max number of intervals to extract */ 
#define MAX_NUM_STR          10   /* max length of an integer string    */
#define MAX_FILENAME_LEN    255   /* max length of unix filename        */
#define MAX_LINE_LEN       1024   /* max length of input lines          */
#define MAX_DESC_LEN       1024   /* max length of optional desc.       */
#define NIL                ((void *) 0)


typedef struct pair_n {  int   a;   
                         int   b;
                         char *desc;
                      } PAIR;

PAIR intervals[MAX_NUM_INTS+1];   /* user-specified intervals go into this */
int  num_ints = 0;                /* table                                 */ 

char  header[MAX_LINE_LEN+1];     /* fasta header is a global */

void  usage()
   {
    fprintf( stderr, "\n\
Usage:       intervals -i n-m,[n-m,...]   [sequence-file]                 \n\
             intervals -f file  [sequence-file]                           \n\
                                                                          \n\
             There can be only one input sequence, and it must be in      \n\
             fasto format.  If no file is specified, or is \"-\", stdin is\n\
             read.                                                        \n\
                                                                          \n\
             One of the two mutualy exclusive options -i or -q must be    \n\
             used to specify the desired intervals.                       \n\
                                                                          \n\
             The first sequence position is 0.  Intervals n-m start and   \n\
             include n and stop at m, but do not include m.  They may     \n\
             overlap.                                                     \n\
                                                                          \n\
Options:        -i n-m,...  comma-seperated intervals are specified on the\n\
                         command line.  Each inteval is a n-m pair of     \n\
                         hyphen-separated integers                        \n\
                                                                          \n\
                -f file    read intervals from file, each pair on one line\n\
                           (- or space separated).                        \n\
                                                                          \n\
                                                                          \n\
" );
    exit( 1 );
   }


void  *malloc_safely( size_t n_bytes )         /* malloc(), or die trying */

   {
    void *p;

    if ( !(p = malloc( n_bytes )) )
       {
        fprintf( stderr, "failed to malloc %d bytes\n", (int) n_bytes );
        exit( errno );
       }
    return( p );
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

/* close_file() closes the file unless its stdin */

void  close_file( FILE *f )
   {
    if ( f != stdin )
        fclose( f );
   }

/* convert numeric string between pointers s and e into an int and */
/* return it, or die trying                                        */

int  get_num( char *s, char *e ) 
   {
    static char nums[MAX_NUM_STR+1];
    int    l;
    l = (int)(e - s);
    if ( l > MAX_NUM_STR )
       {
        fprintf( stderr, "interval position has too many digits\n" );
        usage();
       }
    strncpy( nums, s, l );
    nums[l] = '\0';
    return( (int) strtol( nums, NULL, 10 ) );
   }


void  enter_interval( int x, int y, char *desc )  /* enter interval into global table */
   {
    if ( num_ints < MAX_NUM_INTS )
       {
        if ( x < y )
           {  
            intervals[num_ints].a = x;
            intervals[num_ints].b = y;
           }
        else
           {  
            intervals[num_ints].a = y;   /* force low, high order */
            intervals[num_ints].b = x;
           }
        intervals[num_ints].desc = strndup( desc, MAX_DESC_LEN );
        num_ints++; 
       }
    else
       {
        fprintf( stderr, "max number of intervals exceeded: %d\n",
                          MAX_NUM_INTS );
        exit( 1 );
       }
   }

#ifdef DEBUG

void  print_intervals( void )
   {
    int  i;
    printf( "%d intervals:\n", num_ints );
    for ( i = 0; i <= num_ints; i++ )
        printf( "%d: %d %d\n", i, intervals[i].a, intervals[i].b );
   }

#endif

/* parse string of integer pairs a-b,u-v,n-m,... and fill global */ 

#define ST_NEW    0
#define ST_NUM1   1
#define ST_HYPHEN 2
#define ST_NUM2   3

void  get_intervals_from_string( char *s )
   {
    int state = ST_NEW;  /* state machine state variable, takes above values */
    char *np;            /* points to beginning of current number string     */
    int   a;             /* integer value of first number in pair            */
 
    while ( *s )
       {
        if ( *s != ' ' )          /* this just throws out spaces -work out */
           {                      /* a better system later */
            switch ( state )
               {
                case  ST_NEW:     if ( isdigit( *s ) )
                                     {
                                      np = s;
                                      state = ST_NUM1;
                                     }
                                  else
                                      usage();
                                  break;

                case  ST_NUM1:    if ( *s == '-' )
                                     {
                                      a = get_num( np, s );
                                      state = ST_HYPHEN;
                                     } 
                                  else if ( ! isdigit( *s ) )
                                      usage();
                                  break;
                case  ST_HYPHEN:  if ( isdigit( *s ) )
                                     {
                                      np = s;
                                      state = ST_NUM2;
                                     }
                                  else
                                      usage();
                                  break;
                case  ST_NUM2:    if ( *s == ',' )
                                     {
                                      enter_interval( a, get_num( np, s ), "" );
                                      state = ST_NEW;
                                     }
                                  else if ( ! isdigit( *s ) )
                                      usage();
                                  break;
               }
           }
        s++;
       }
    if ( state == ST_NUM2 )
        enter_interval( a, get_num( np, s ), "" );
    else if ( ! isdigit( *s ) )
        usage();
   }

/* This variation reads a-b interval pairs from a file, one pair  */
/* per line.                                                      */

void  read_intervals_from_file ( char *filename )
   {
    FILE *f;
    static char line[MAX_LINE_LEN+1];
    char   *s;
    int    nvals;
    int    a, b;
    static char desc[MAX_DESC_LEN+1];

    f = open_file( filename );
    while ( fgets( line, MAX_LINE_LEN, f ) )
       {
        for ( s = line; *s != '\0' && *s != '\n'; s++ ) /* replace first occur*/
            if ( *s == '-' )                            /* ance of '-'      */
               {
                *s = ' ';
                break;
               }
        nvals = sscanf( line, "%d %d %[^\n]s", &a, &b, desc );
        if ( nvals == 3 )
            enter_interval( a, b, desc );
        else if ( nvals == 2 )
            enter_interval( a, b, "" );
        else if ( nvals != 0 )
           {
            fprintf( stderr, "bad line in intervals file %s: %s\n",
                     filename, line );
            close_file( f );
            exit( 1 );
           }
       }
    close_file( f );
   } 

/* read command line arguments, set globals accordingly */

void  parse_args( int argc, char **argv, char *file )
   {
    extern char *optarg;
    extern int   optind;
    int          c;

    while ( (c = getopt( argc, argv, "i:f:h" ) ) != -1 )
        switch( c )
           {
            case  'i':   get_intervals_from_string( optarg );
                         break;
            case  'f':   read_intervals_from_file( optarg );
                         break;
            case  'h':   usage();
            default:     usage();
           }
    argc -= optind;
    argv += optind;
    if ( argc == 1 )
        strncpy( file, *argv, MAX_FILENAME_LEN );
    else if ( argc == 0 )
        strcpy( file, "-" ); 
    else
        usage();
   }

             


int  pair_cmp( const void *x, const void *y ) /*used for qsort() on intervals*/
   {
    if ( ((PAIR *) x)->a == ((PAIR *) y)->a )
        return( 0 );
    else
        if (  ((PAIR *) x)->a < ((PAIR *) y)->a )
            return( -1 );
        else
            return( 1 );
   }


int  is_sequence_char( int c )   /* returns 1 unlesss c is whitespace */
   {
    return( !( c == ' ' || c == '\n' || c == '\t' || iscntrl( c ) ) );
   }


/* print sequence s in fasta form*/
void  output_fasta( char *s, int a, int b, char *desc )
   {                                        /* (and add ":a-b" to end of hdr)*/
    int         i, n;

    printf( ">%s %d-%d %s\n", desc, a, b, header+1 );
    n = strlen( s );
    for ( i =  50; i <= n; i += 50, s += 50 )
        printf( "%50.50s\n", s );
    if ( *s != '\0' )
        printf( "%s\n", s );
   }

/*************************************************************************/
/* Output Queues structure:  a linked list of nodes of type QUEUE, one   */
/* for each "active" interval - as the pos counter moves through the     */
/* main sequence, it passes through the specified intervals, some which  */
/* may overlap.  open_output_queue() is invoked each time a new interval */
/* is entered, which adds a new queue to the list.  enter_queue() adds   */
/* characters to the individual output queues, and flush_and_close()     */
/* prints out the queue sequences and removes the queues                 */
/*************************************************************************/

typedef struct q_p { 
                     char        *seq; /* points to a sequence buffer */
                     int          a;   /* from corresponding intervals entry */
                     int          b;
                     char        *desc;
                     struct q_p  *next; /* next node in linked list */
                   } QUEUE;

QUEUE *queue_top = NIL;


void  dump_queues( void )  /* for debuging */
   {
    QUEUE *p;
    printf( "queues now:\n" );
    for ( p = queue_top;  p != NIL;  p = p->next )
        printf( "  %d-%d %s\n", p->a, p->b, p->desc ); 
   }

void  open_output_queue( int i )
   {
    QUEUE *q;

    q = malloc_safely( sizeof( QUEUE ) );
    q->a = intervals[i].a;                        /* since we already know  */
    q->b = intervals[i].b;                        /* the string size, we can*/
    q->desc = intervals[i].desc;
    q->seq = malloc_safely( q->b - q->a + 1 );    /* allocate string storage*/
    q->next = queue_top;                          /* finally, add to the top*/
    queue_top = q;                                /* of the list            */
   }

int queues_still_open( void )                     /* return 1 if there are */
   {                                              /* entries in the list   */
    return( queue_top != NIL ); 
   } 

int  past_interval( int pos, QUEUE *q )           /* returns 1 if we've have*/
   {                                              /* exited the interval    */
    return( pos >= q->b );
   }

                                             
void  enter_queue( QUEUE *q, int pos, int c )     /* add c to sequence buffer*/
   {                                              /* for interval q          */
    q->seq[pos - q->a] = c;
   }

/* print out the sequence, for interval q, remove from list, and deallocate */
/* all storage.  As a convenience, this also returns a pointer to the next  */
/* node, because q->next is meaningless once q has been free'd              */

QUEUE *flush_and_close( QUEUE *q, int pos )
   {
    QUEUE *p;                                    

    if ( q == NIL )
       { 
        fprintf( stderr, "Internal error: flush_and_close() received NIL\n");
        exit( 1 );
       }
    enter_queue( q, pos, '\0' );                 /* need to terminate string */
    output_fasta( q->seq, q->a, q->b, q->desc ); /* and then print it out and*/
    free( q->seq );                              /* then free up its memory  */
    if ( q == queue_top )                        /* if q is at the top, then */
        p = queue_top = q->next;                 /* this is what we do,      */
    else                                         /* otherwise we have to     */
       {                                         /* search the list for it   */
        for ( p = queue_top; p->next != NIL && p->next != q; p = p->next ) 
            ;
        p->next = q->next;                       /* and make a shunt around  */
        p = q->next;                             /* around it, and then */ 
       }
    free( q );                                   /* free it up, and also     */
    return( p );                                 /* return a pointer to next */
   } 



                                   /****************/
                                   /* Main Program */
                                   /****************/


int  main( int argc, char **argv )
   {
    int          c;
    int          pos;
    static char  filename[MAX_FILENAME_LEN+1];
    FILE        *f;
    int          next_int = 0;
    QUEUE       *q;

    parse_args( argc, argv, filename );  /* this fills the interval table */

    qsort( intervals, num_ints, sizeof( PAIR ), pair_cmp );  /* sort it */

#ifdef DEBUG
    print_intervals();
#endif   
    f = open_file( filename );                      /* we assume the file */ 
    fgets( header, MAX_LINE_LEN, f );               /* has only one fasta */ 
    header[strlen(header)-1] = '\0';                /* format sequence    */

    pos = -1;                                       /* keep going while there*/
    while ( (queues_still_open() || next_int < num_ints) /* are open queues, */
             && (c = fgetc( f )) )                       /* more intevls, and*/
                                                         /* and more input   */
        if ( is_sequence_char( c ) )                     /* ignore whitespace*/
           {
            pos++;                                       /* position counter */

            /* have we now entered any new intervals? if so, open queues */

            while ( (next_int < num_ints) && pos >= intervals[next_int].a )
                open_output_queue( next_int++ );

            q = queue_top;                          /* each open queue, if  */
            while ( q != NIL )                      /* we've exited, print  */
                if ( past_interval( pos, q ) )      /* it and delete it     */
                    q = flush_and_close( q, pos );
                else                                /* otherwise append this*/
                   {                                /* character to the seq */ 
                    enter_queue( q, pos, c );       /* sequence buffer and  */
                    q = q->next;                    /* proceed to the next  */
                   }                                /* queue                */
           }
   
    close_file( f );
   }


