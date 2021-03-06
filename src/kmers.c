/* Program:     kmers.c                                                      */
/* Programmer:  Sean R. McCorkle                                             */
/*              Biology Dept. Brookhaven National Laboratory                 */
/* Language:    C                                                            */
/*                                                                           */
/* Description: counts k-mer frequencies in DNA sequences                    */
/*                                                                           */
/* Notes:       Histogram indices are constructed from packed two-bit        */
/*              encodings of the nucleotides (A->0, C->1, G->2, T->3) so     */
/*              maximum value of k, MAX_K, is limited to 16 for 32 bit words */
/*                                                                           */
/*              A further limitation, perhaps more severe, is the size of the*/
/*              histogram, which is 4^k.  (see MAX_HIST_LEN)                 */
/*                                                                           */
/*              This is case insensitive. A=a, C=c, G=g, T=t at all times.   */
/*                                                                           */
/*              Ignores k-mers with any ambiguity codes or any other         */
/*              characters which are not nucleotides A, C, G, T.             */
/*                                                                           */
/*                                                                           */
/* Compiling:   cc -O -o kmers kmers.c should do the job.                    */
/*                 (or cc -O3 if you prefer)                                 */
/*                                                                           */
/*****************************************************************************/

#include <errno.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <unistd.h>

static char rcsvers[] = "$Revision: 1.3 $";

#define MAX_K              10  /* this is constrained by the limit on 4^k */
#define MAX_HIST_LEN  1048576  /* must be 4^MAX_K */
#define N_ASCII           128  /* size of ascii arrays (7 bits) */
                               /* hmm, not sure what will happen with 8 bits*/


/* globals set by command line options */

int  word_size          = 2;  /* the k in k-mers; set by -k  */
int  report_by_sequence = 0;  /* set by -s option */
int  print_total        = 0;  /* set by -T */
int  verbose            = 0;  /* set by -v */

#define  A          65       /* ASCII codes for nucleotides */
#define  C          67
#define  G          71
#define  T          84
#define  NIL        '_'      /* initially empty state */
#define  LOWER_OFF  32       /* offset for lowercase version of letters */

unsigned int  real_nts[] = { A, C, G, T, 0 };
                             /* this order ensures bitwise complimentarity */
                             /* which int_comp() uses */


#define  MAX_HEADER_LEN    16384          /* was 1024, needed to handle Trinity*/
char     header[MAX_HEADER_LEN+1];        /* holds a fasta header */

/* these arrays act as character functions - each maps an ascii value */
/* to a function value                                                */

unsigned int  twobit[N_ASCII];            /* maps nucleotide char to 2 bits */
unsigned int  no_count[N_ASCII];          /* 1 => char is not to be counted */
unsigned int  is_not_space[N_ASCII];      /* 1 if whitespace */

/* these arrays map k-mer (in compressed integer form) to a value. */

unsigned long int  counts[MAX_HIST_LEN];
char              *kmer_seq[MAX_HIST_LEN];
unsigned long int  rc_map[MAX_HIST_LEN];  /* maps index to index of rc */

/* more globals */

int                n_kmers;            /* determined by k:  4^k              */
int                n_seqs = 0;         /* number of sequences processed      */

/* these next 5 form the circular buffer apparatus, which is reset           */
/* by reset_cbuff() and updated by insert_cbuff()                            */

int                n_nocounts;         /* tracks number of non-nucleotide    */
                                       /* chars currently in circular buffer */
unsigned long int  w_mask;             /* bit mask for controlling circular  */
                                       /* index                              */
unsigned long int  cbuff_w;            /* circular buffer index; contains    */
                                       /* k-mer in a compressed two-bit form,*/
                                       /* i.e. 0x1e corresponds to 3-mer CTG */
char               cbuff[MAX_K];       /* ASCII char circular buffer, for    */
                                       /* tracking outgoing characters       */
int                next_cbuff = 0;     /* next position in cbuff to accept   */


/* version() - print program name and version */

void  version( void )
   {
    char *v;
    for ( v = rcsvers; *v && *v != ' '; v++ )
        ;
    v++;
    v[strlen(v)-1] = '\0';
    fprintf( stderr, "\nkmers        version %s\n", v );
   }


/* usage() - print help */

void usage( void )
   {
    fprintf( stderr, " \n\
Usage:       kmers [-k<n>] [-hsTvV]  [seq-file ... ]                      \n\
                                                                          \n\
             where [seq-files] are in FASTA format.  The name \"-\" means \n\
             stdin.  stdin is scanned it no filemames are specified.      \n\
                                                                          \n\
Options:     -k<n>   count k-mers of size <n>                             \n\
             -s      print counts for each sequence                       \n\
             -T      print a total of dimer counts (1-direction)          \n\
             -v      verbose mode                                         \n\
             -V      print version                                        \n\
             -h      print help message                                   \n\
                                                                          \n\
" );

   }


/* extract command line arguments and options */

void parse_args ( int argc, char **argv, int *nfiles, char ***filenames )
   {
    static char *stand_in[] = { "-", (char *) 0 };
    int    c;
    char  *endptr;

    while ( (c = getopt( argc, argv, "k:TsvVh")) != -1 )
        switch ( c )
           {
            case 'k':  word_size = strtol( optarg, &endptr, 10 );
                       if ( errno == 0 )
                          {
                           if ( endptr == optarg )
                               usage();
                           else if ( word_size<1 || word_size>MAX_K )
                              {
                               fprintf( stderr, 
                                          "k must be in range 1-%d\n", MAX_K);
                               exit( 1 );
                              }
                          }
                       else
                          { 
                           perror( "k argument conversion error" );
                           exit( errno ); 
                          }
                       break;
            case 's':  report_by_sequence = 1;  break;
            case 'T':  print_total = 1;         break;
            case 'v':  verbose = 1;             break;
            case 'V':  version();               exit(0);
            case 'h':  version();
                       usage();
                       exit(0);
            default:   usage();                 exit(1);
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
        *filenames = stand_in;
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
            exit( errno );
           }
   }


/*************************************************/
/* close_file() closes the file unless its stdin */
/*************************************************/

void  close_file( FILE *f )
   {
    if ( f != stdin )
        fclose( f );
   }



#ifdef DEBUG
void  dumps( char *msg, unsigned int a[] )
   {
    int i;
    printf( "%s\n", msg );
    for ( i = 0; i < N_ASCII; i++ )
        if ( a[i] )
            printf( "  %3d %c %d\n", i, i, a[i] );
    }
#endif


/* this returns as an ASCII string, the sequence corresponding to     */
/* the given compressed 2-bit index.  This doesn't allocate any new   */
/* memory - that must be handled in the calling routine               */

char  *int2seq( int n )
   {
    static char seq[MAX_K+1];
    int         j;

    seq[word_size] = '\0';
    for ( j = word_size-1;  j >= 0;  j-- )
       {
        seq[j] = real_nts[n & 0x3];
        n >>= 2;
       }
    return( seq );
   }


/* this returns the compressed bit index of the k-mer which is the */
/* reverse complement of the given k-mer index                     */

unsigned int  intcomp ( int n )
   {
    int          j;
    unsigned int comp = 0x0;

    for ( j = 0; j < word_size; j++ )
       {
        comp = (comp << 2) | ( n & 0x3 );
        n >>= 2;
       }
    /* at this point, comp corresponds to the reverse sequence, so */
    /* simply bitwise-complement and return */

    return( w_mask & ~comp );
   }


/* init() prefills ALL the character function arrays and also         */
/* prepares the kmer-index to sequence map and reverse complement map */

void  init( void )
   {
    int i;

    /* fill twobits[] */
    bzero( twobit, N_ASCII * sizeof( unsigned int ) );
    for ( i = 0; i < 4; i++ )
         twobit[real_nts[i]] = twobit[real_nts[i]+LOWER_OFF] = i;
#ifdef DEBUG
    dumps( "twobit", twobit );
#endif

    /* fill no_count[] - don't count anything ...*/
    for ( i = 0; i < N_ASCII; i++ )
        no_count[i] = 1;
    /*... except real nucleotides */
    for ( i = 0; i < 4; i++ )
         no_count[real_nts[i]] = no_count[real_nts[i]+LOWER_OFF] = 0;
#ifdef DEBUG
    dumps( "no_count", no_count );
#endif

    /* is_not_space[c] is true for any chars which are not control  */
    /* chars or whitespace (make whitespace explicit here!)         */

    for ( i = 0; i < N_ASCII; i++ )
        if ( i >= 33 && i <= 126 && i != ' ' && i != '\t' && i != '\n' )
            is_not_space[i] = 1;
        else
            is_not_space[i] = 0;
#ifdef DEBUG
    dumps( "is_not_space", is_not_space );
#endif

    w_mask = 0;
    n_kmers = 1;
    for ( i = 0; i < word_size; i++ )
       {
        w_mask = (w_mask << 2) | 0x3;
        n_kmers *= 4;
       }

#ifdef DEBUG
    printf( "word_size is %d, w_mask is %x, n_kmers is %d (0x%x)\n", 
             word_size, w_mask, n_kmers, n_kmers );
#endif

    for ( i = 0; i < n_kmers; i++ )
       {
        kmer_seq[i] = strdup( int2seq( i ) );
#ifdef DEBUG
        printf( "%10d 0x%08x [%s]\n", i, i, kmer_seq[i] );
#endif
       }    
    for ( i = 0; i < n_kmers; i++ )
       {
        rc_map[i] = intcomp( i );
#ifdef DEBUG
        printf( "%10d 0x%08x [%s] -> %10d 0x%08x [%s]\n", 
            i, i, kmer_seq[i], rc_map[i], rc_map[i], kmer_seq[rc_map[i]] ); 
#endif
       }    
   }


/* clear k-mer count histogram */

void  clear_counts( void )
   {
    bzero( counts, n_kmers * sizeof( unsigned long int ));
   }


#ifdef DEBUG
void  dump_cbuff( void )
   {
    int i;
    printf( "cbuff_w %x  n_nocounts = %d\n", cbuff_w, n_nocounts );
    for ( i = 0; i < word_size; i++ )
        printf( "   %d: %c%s\n", i, cbuff[i], (i == next_cbuff)? " *" : "" );
   }
#endif


void  reset_cbuff( void )
   {
    int  i;

#ifdef DEBUG
    printf( "reset cbuff\n" );
#endif
    next_cbuff = 0;
    for ( i = 0; i < MAX_K; i++ )
       cbuff[i] = NIL;
    cbuff_w = w_mask & 0xffffffff;
    n_nocounts = word_size;   /* cbuff starts out filled with NILs */
   }


/* this clears and resets everything */

void  reset( void )
   {
    reset_cbuff();
    clear_counts();
#ifdef DEBUG
    dump_cbuff();
#endif
   }


/* insert one character (non-whitespace) into the cbuff */

void  enter_cbuff( char c )
   {
    if ( no_count[ cbuff[next_cbuff] ] ) /* is outgoing char ambiguity code? */
        n_nocounts--;                    /* if so, reduce count */
    if ( no_count[c] )                   /* is incoming char ambiguity code? */
        n_nocounts++;                    /* if so, increment count */
    cbuff[ next_cbuff++ ] = c;
    if ( next_cbuff >= word_size )
        next_cbuff = 0;
    cbuff_w = w_mask & ( (cbuff_w << 2) | twobit[c] );
   }


/* use the current circular buffer compressed index and increment */
/* the counts for that index.  (note -this is only doing the top  */
/* strand - we'll infer the bottom strand counts in report()      */

void incr_counts( void )
   {
    if ( n_nocounts == 0 )
       { 
#ifdef DEBUG
        printf( "incr counts[%d] (%x)\n", cbuff_w, cbuff_w ); 
#endif
        counts[cbuff_w]++; 
       }
    else if ( n_nocounts < 0 )
       {
        fprintf( stderr, "Yeow negative n_nocounts: %d\n", n_nocounts );
        exit( 1 );
       }
   }

/* calculate total nt by adding up appropriate entries from array    */
/* counts[].  For now, we will not include any ambiguities unless -a */
/* set */

unsigned long int  compute_total( void )
   {
    int  i;
    unsigned long int  total = 0;

    for ( i = 0; i < n_kmers; i++ )
        total += counts[i];
    return( total );
   }

/* here's the report.  HEre, we calculate the total (for percentages) */
/* AND where we handle the bottom strand counts (making use of the    */
/* rc_map[] array to guide us to each k-mer's reverse complement      */

void  report( void )
   {
    int i;
    unsigned long int  total;

    if ( report_by_sequence )
        printf( ">%s\n", header );
    total = compute_total();
    for ( i = 0; i < n_kmers; i++ )
       {
        printf( "%s/%s %10d %10d %10d  %10.6f %10.6f %10.6f\n",
                 kmer_seq[i], kmer_seq[rc_map[i]], 
                 (int) counts[i], (int) counts[rc_map[i]], (int) (counts[i] + counts[rc_map[i]]),
                 100.0 * ( (double) counts[i] / (double) total ),
                 100.0 * ( (double) counts[rc_map[i]]  / (double) total ),
                 50.0 * ( (double) counts[i] + (double) counts[rc_map[i]] )
                         / (double) total
                );
       }    
    if ( print_total )
        printf( "total: %12.0f\n", (double) total );
   }

/* This routine is reponsible for processing one file, which has already */
/* been opened for it  */

void  count_kmers( FILE *f )
   {
    int c;
    int last_c = '\n';
    int i;

    while ( (c = fgetc( f )) != EOF )
       {
        if ( last_c == '\n' && c == '>' )  /* is this a new sequence header? */
           {                               /* if so, then if we're reporting */
            if ( n_seqs++ > 0 && report_by_sequence  )  /* each sequence,then*/
                {                                       /* report and clear */
                 report();                              /* counters */
                 reset();
                }
            else                          /* in any case, we need to reset */
                reset_cbuff();            /* the circular buffer apparatus */
 
            i = 0;                        /* and grab the header in any case */
            while ( i < MAX_HEADER_LEN &&(c = fgetc( f )) != EOF && c != '\n' )
                header[i++] = c;
            if ( i >= MAX_HEADER_LEN )
               {
                header[MAX_HEADER_LEN] = '\0';
                fprintf( stderr, "header too long: %s\n", header );
                exit( 1 );
               }
            header[i] = '\0';
            if ( verbose )
                fprintf( stderr, ">%s\n", header );
           }
        else if ( is_not_space[c] )      /* otherwise, if not  whitespace  */
           {                             /* then process, by inserting in  */
            enter_cbuff( c );            /* into the circular buffer, then */
            incr_counts();               /* updating the counters          */
           }

        last_c = c;                      /* don't forget this! */
       }
   }

                                /****************/
                                /* Main Program */
                                /****************/



int main ( int argc, char **argv )
   {
    int           nfiles;
    int           i;
    static char **filenames;
    FILE         *f;

    parse_args( argc, argv, &nfiles, &filenames );
    init();
    reset();

    for ( i = 0; i < nfiles; i++ )     /* for each file */
       { 
        if ( verbose )
            fprintf( stderr, "file: %s\n", filenames[i] );
        f = open_file( filenames[i] );
        count_kmers( f );
        close_file( f );
       }
    report();
   }

