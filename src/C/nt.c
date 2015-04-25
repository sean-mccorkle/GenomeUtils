/* Program:     nt.c                                                         */
/* Programmer:  Sean R. McCorkle                                             */
/*              Biology Dept. Brookhaven National Laboratory                 */
/* Language:    C                                                            */
/*                                                                           */
/* Description: counts single nucleotide frequencies in DNA sequences        */
/*                                                                           */
/*                                                                           */
/* Compiling:   cc -O -o nt nt.c should do the job.                          */
/*                 (or cc -O3 if you prefer)                                 */
/*                                                                           */
/* $Id: nt.c,v 1.5 2012/09/27 21:06:21 seanmccorkle Exp seanmccorkle $           */
/*****************************************************************************/

#include <errno.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

static char nt_rcs_id[] =
    "$Id: nt.c,v 1.5 2012/09/27 21:06:21 seanmccorkle Exp seanmccorkle $";
static char rcsvers[] = "$Revision: 1.5 $";

int  show_ambs          = 0;  /* set to 1 by -a option;  to 2 by -A */
int  case_sensitive     = 0;  /* set by -c option */
int  gc_content         = 0;  /* set by -g option */
int  suppress_zeros     = 0;  /* set by -z option */
int  report_by_sequence = 0;  /* set by -s option */
int  show_nts           = 1;  /* cleared by -n */
int  show_counts        = 1;  /* cleared by -P */
int  show_percents      = 0;  /* set by -p, or -P */
int  show_total         = 0;  /* set by -T */

#define  MAX_HEADER_LEN   16384  /* was 1024, but encountered a big one once */
char     header[MAX_HEADER_LEN+1];

int  n_seqs = 0;


/* version() - print program name and version */

void  version( void )
   {
    char *v;
    for ( v = rcsvers; *v && *v != ' '; v++ )
        ;
    v++;
    v[strlen(v)-1] = '\0';
    fprintf( stderr, "\nnt           version %s\n", v );
   }


/* usage() - print help */

void usage( void )
   {
    fprintf( stderr, " \n\
Usage:       nt  [-aAcgnhpPsTVz]  [seq-file ... ]                         \n\
                                                                          \n\
             where [seq-files] are in FASTA format.  The name \"-\" means \n\
             stdin.  stdin is scanned it no filemames are specified.      \n\
                                                                          \n\
Options:     -a      treat all ambiguities as \"N\", count & display      \n\
             -A      count and display all ambiguities (N,M,R,Y,etc)      \n\
             -c      case sensitive: count upper & lower case separately  \n\
             -g      count and display G/C and A/T freqs. (case ins.)     \n\
             -n      don't print nucleotides, just counts/percents        \n\
             -p      print percentages along with total counts            \n\
             -P      print only percentags (not total counts)             \n\
             -s      print counts/percents for each fasta sequence        \n\
             -T      print total nucleotide count                         \n\
             -z      suppress counts zero in output                       \n\
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

    while ( (c = getopt( argc, argv, "aAcghnpPsTVz")) != -1 )
        switch ( c )
           {
            case 'a':  show_ambs = 1;           break;
            case 'A':  show_ambs = 2;           break;
            case 'c':  case_sensitive = 1;      break;
            case 'g':  gc_content = 1;          break;
            case 'n':  show_nts = 0;            break;
            case 'p':  show_percents = 1;       break;
            case 'P':  show_percents = 1;       
                       show_counts = 0;         break;
            case 's':  report_by_sequence = 1;  break;
            case 'T':  show_total = 1;          break;
            case 'z':  suppress_zeros = 1;      break;
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


/* close_file() closes the file unless its stdin */

void  close_file( FILE *f )
   {
    if ( f != stdin )
        fclose( f );
   }


unsigned long int  counts[256];
unsigned long int  total_nt;

void  clear_counts( void )
   {
    bzero( counts, 256 * sizeof( unsigned long int ));
   }



/* ascii codes for convenience */

#define  NON_PRINT  32
#define  A          65
#define  C          67
#define  G          71
#define  T          84
#define  N          78   /* amb code A, C, G, or T */
#define  M          77   /* amb code A or C */
#define  R          82   /*          A or G */
#define  W          87   /*          A or T */
#define  S          83   /*          C or G */
#define  Y          89   /*          C or T */
#define  K          75   /*          G or T */
#define  V          86   /*          A, C or G */
#define  H          72   /*          A, C or T */
#define  D          68   /*          A, G or T */
#define  B          66   /*          C, G or T */
#define  X          88
#define  LOWER_OFF  32  /* offset for lowercase version of letters */

#define  GC        0xfc  /* internal codes for output  for -g option */
#define  AT        0xfd

int  real_nts[] = { A, C, G, T, 0 };
int  ambiguities[] = { N, M, R, W, S, Y, K, V, H, D, B, X, 0 };  
                     /* Not sure how to handle X actually.  ambiguity? */


/* calculate total nt by adding up appropriate entries from array    */
/* counts[].  For now, we will not include any ambiguities unless -a */
/* set.  Result goes into global total_nt                            */

void  compute_total( void )
   {
    int  k;
   
    total_nt = 0;
    for ( k = 0; real_nts[k] > 0; k++ )
        total_nt += counts[ real_nts[k] ] + counts[ real_nts[k] + LOWER_OFF ];
    if ( show_ambs > 0 )
        for ( k = 0; ambiguities[k] > 0; k++ )
            total_nt += counts[ ambiguities[k] ] + 
                                    counts[ ambiguities[k] + LOWER_OFF ];
   }


/* this prints the counts and or percentags for one nucleotide, */
/* for any of the output format options                         */

void  print_count( int nt, int cnt )
   {
    double percnt;

    if ( show_nts )
       {
        if ( nt == GC )
            printf( "G/C: " );
        else if ( nt == AT )
            printf( "A/T: " );
        else
            printf( "%c: ", nt );
       }
    if ( show_counts )
        printf( "%10d ", cnt );
    if ( show_percents )
       {
        if ( total_nt > 0 )
            percnt = (100.0 * cnt) / total_nt;
        else
            percnt = -1.0;
        printf( "%6.2lf ", percnt );
       }
    if ( ! report_by_sequence )
        printf( "\n" );
   }


/* handle the reporting of one nucleotide, for any output format option */

void  report_nt( int nt )
   {
    int  i;
    int  tot;
    int  tot_low;

    if ( nt == GC )
        print_count( nt, counts[G] + counts[G+LOWER_OFF] + 
                         counts[C] + counts[C+LOWER_OFF] );
    else if ( nt == AT )
        print_count( nt, counts[A] + counts[A+LOWER_OFF] + 
                         counts[T] + counts[T+LOWER_OFF] );
    else if ( nt == N && show_ambs == 1 )
       {   /* count all ambs as N, and respect the case sensitivity */
           /* hopefully these are low counts, i.e. no overflow */
        tot = 0;
        if ( case_sensitive )
           {
            tot_low = 0;
            for ( i = 0; ambiguities[i]; i++ )
               {
                tot += counts[ ambiguities[i] ];
                tot_low += counts[ ambiguities[i] + LOWER_OFF ];
               }
            print_count( nt, tot );
            print_count( nt+LOWER_OFF, tot_low );
           }
        else
           {
            for ( i = 0; ambiguities[i]; i++ )
                tot += counts[ ambiguities[i] ] + 
                         counts[ ambiguities[i] + LOWER_OFF ];
            print_count( nt, tot );
           }
       }
    else if ( case_sensitive )
       {
        print_count( nt, counts[nt] );
        print_count( nt+LOWER_OFF, counts[nt+LOWER_OFF] );
       }
    else
        print_count( nt, counts[nt] + counts[nt+LOWER_OFF] );
   }


/* output results, for any output format option */

void  report( void )
   {
    int i;

    compute_total();
    if ( gc_content )
       {
        report_nt( GC );
        report_nt( AT );
       }
    else
       {
        report_nt( A );
        report_nt( C );
        report_nt( G );
        report_nt( T );
        if ( show_ambs == 1 )
            report_nt( N );
        else if ( show_ambs > 1 )
            for ( i = 0; ambiguities[i]; i++ )
                report_nt( ambiguities[i] );
            
      }
    if ( show_total )
        printf( "total:  %ld\n", total_nt );
    if ( report_by_sequence )
        printf( "%s\n", header );
   }

/* count the nucleotides in one file, which is already opened     */
/* If report_by_sequence (-s option), output counts and reset for */
/* each DNA sequence in the file                                  */

void  count_nts( FILE *f )
   {
    int c;
    int last_c = '\n';  /* last char helps identify fasta header record*/
    int i;

    while ( (c = fgetc( f )) != EOF )
       {
        if ( last_c == '\n' && c == '>' )
           {
            if ( n_seqs++ > 0 && report_by_sequence  )
                {
                 report();
                 clear_counts();
                }
 
            i = 0;
            while ( i < MAX_HEADER_LEN &&(c = fgetc( f )) != EOF && c != '\n' )
                header[i++] = c;
            if ( i >= MAX_HEADER_LEN )
               {
                header[MAX_HEADER_LEN] = '\0';
                fprintf( stderr, "header too long: %s\n", header );
                exit( 1 );
               }
            else
               header[i] = '\0';
           }
        else
            counts[c & 0x00ff]++;

        last_c = c;
       }
   }

                             /****************/
                             /* Main Program */
                             /****************/


main ( int argc, char **argv )
   {
    int           nfiles;
    int           i;
    static char **filenames;
    FILE         *f;

    parse_args( argc, argv, &nfiles, &filenames );
    clear_counts();
    for ( i = 0; i < nfiles; i++ )
       { 
        f = open_file( filenames[i] );
        count_nts( f );
        close_file( f );
       }
    report();
   }

