/* File:         seqlib.h                                                 */
/* Programmer:   Sean R. McCorkle                                         */
/* Lanuage:      C                                                        */
/* Description:  Definitions needed for seqlib routines.                  */
/*               query routines.                                          */
/*                                                                        */
/* $Id: seqlib.h,v 3.5 2007/11/13 19:29:56 mccorkle Exp mccorkle $        */
/**************************************************************************/

#include <stdlib.h>

                         /*********************/
                         /* General Constants */
                         /*********************/

#define  STD_ALLOWED      "acgtmrwsykvhdbnxACGTMRWSYKVHDBNX"
#define  STD_IGNORED      " \t\n"

#define  MAX_HDR_LENGTH     1024   /* max length of fasta >header line       */
#define  MAX_SEQ_LENGTH    12000   /* max length of sequence we can tolerate */
#define  MAX_LINE_LENGTH     256   /* max length for input lines             */
#define  MAX_HEADER_LINES    256   /* max # of seq file hdr & comment lines  */

#define  ERROR                -2   /* load_seq() return codes                */
#define  NUM_CHARS           256   /* size of char set lookup tables         */


#define  X                    -1   /* not an allowed character               */

                         /*******************/
                         /* Data Structures */
                         /*******************/

struct align_stats_s {
                       int  top_beg_off;       /* leading inserts in seq1   */
                       int  top_end_off;       /* trailing inserts in seq1  */
                       int  bot_beg_off;       /* leading inserts in seq2   */
                       int  bot_end_off;       /* trailing inserts in seq2  */
                       int  num_ok;            /* number of exact matches   */
                       int  num_wrong;         /* number of subst. errors   */
                       int  num_indels;        /* number of indels          */
                       int  num_seq1_inserts;  /* number of inserts in seq1 */
                       int  num_seq2_inserts;  /* number of inserts in seq2 */
                       int  num_ambs;          /* number of ambiguity substs*/
                       int  num_seq1_ambs;     /* number where ambs in seq1 */
                       int  num_seq2_ambs;     /* number where ambs in seq2 */

                       /* these are counts for those which occur in aa  */
                       /* translation changes                           */

                       int  num_wrong_x;        /* number of subst. errors   */
                       int  num_indels_x;       /* number of indels          */
                       int  num_seq1_inserts_x;/* number of inserts in seq1 */
                       int  num_seq2_inserts_x;/* number of inserts in seq2 */

                     };

typedef struct align_stats_s ALIGN_STATS;

                              /*************/
                              /* Externals */
                              /*************/

extern char *progname;   /* declared in sequtils.c - used for error messages */
extern char  seq_map[];  /* declared in sequtils.c - used by SEQ()           */
extern char  agree[16][16];    /* declared in sequtils.c - used by AGREE()   */

                                /**********/
                                /* Macros */
                                /**********/

#define  SEQUENT(x)   ( seq_map[ (x) ] )   /* no check on x - [000 - 0177 ] */
#define  CHARAC(x)    ( (x) < 16 && (x) >= 0 ? ".acgtmrwsykvhdbn"[(x)] : '.' )
#define  AGREE(x,y)   ( agree[x][y] )
#define  ISAMB(x)     ( (x) > 4 )


     /* These are in io.c */

int   open_seqf( char *filename, char *allowed, char *ignored );
int   close_seqf( int s );
int   read_seq( int s, int max_hdr_size, char *hdr,
                int max_seq_size, char *seq);

void  print_seq( FILE *output, char *seq, int start, int stop, int indent,
                 int linelen, int spaces, int xtralines
               );
int   count_char_occurrences( char *set, char *seq, int start, int stop );

FILE *fopen_safely( char *name, char *mode );
void *calloc_safely( size_t nmemb, size_t size, char *errmsg );
void *malloc_safely( size_t size, char *errmsg );


     /* These are in sequtils.c */
void  char_to_sequential( char *characters, char *sequential );
void  sequential_to_char( char *sequential, char *characters );
void  dna_reverse_comp( char *s, char *r );
void  seq_reverse_comp( char *s, char *r );
int   blank( char *s );
char *remove_path( char *path );
char *copystr( int n, char *s );
int   eqchar( char a, char b );
int   eqdna( char a, char b );
void  translate( char *dna, char *aa );

#ifdef METROWERKS
char *strdup( char *s);
#endif

     /* These are in lpa_align.c */
void  lpa_align( char *seq1, int m, char *seq2, int n,
                 int max_align_size, int handle_translation,
                 char *top, char *mid, char *bot, 
                 char *trans_top, char *trans_bot, int *l);
                  
void  print_alignment( FILE *f, char *top_name, char *bot_name,
                       int handle_translation,
                       char *top, char *mid, char *bot, 
                       char *trans_top, char *trans_bot, int k );
void  get_align_stats( char *top, char *mid, char *bot, int l, 
                       int m, int n, ALIGN_STATS *a, int handle_translation, 
                       char *trans_top, char *trans_bot );
void  print_align_stats( ALIGN_STATS *a, int handle_translation );
void  print_diffs( ALIGN_STATS *a, int handle_translation,
                   char *top, char *mid, char *bot, 
                   char *trans_top, char *trans_bot, int n );



