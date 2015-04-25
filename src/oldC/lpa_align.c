/* Module:       lpa_align.c                                              */
/* Programmer:   Sean R. McCorkle                                         */
/* Lanuage:      C                                                        */
/* Description:  least-path string alignment routine and associated stuff */
/*                                                                        */
/* $Id: lpa_align.c,v 3.4 2000/03/22 01:02:50 mccorkle Exp mccorkle $ */
/**************************************************************************/

static char lpa_align_rcs_id[] =
      "$Id: lpa_align.c,v 3.4 2000/03/22 01:02:50 mccorkle Exp mccorkle $";

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "seqlib.h"


                         /************************/
                         /* Macros and Constants */
                         /************************/

typedef  unsigned char  BYTE;

#define  INDEL          5  /* insertion/deletion penalty */
#define  SUB            4  /* substitution error penalty */
#define  AMB            0  /* ambiguity subsititution penalty */

#define  UP             0x04     /* these are the backpointer directions    */
#define  DIAG           0x02     /* they are bit masks so that they can be  */
#define  LEFT           0x01     /* or-ed together to indicate combinations */
#define  UNSEEN         0x7FFFFFFF /* must be positive!!! */

#define  INDEL_IND      '-'      /* used in the mid[] array */
#define  AMB_IND        '|'
#define  WRONG_IND      '*'
#define  OK_IND         ' '
#define  INDEL_SPACER   ' '     /* used in top[] and bot[] array's for indels*/



/* Penalty matrix */

char penalty[16][16] = {
/*          a   c   g   t   m   r   w   s   y   k   v   h   d   b   n */
       0,   X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,  X,
/*a*/  X,   0,SUB,SUB,SUB,AMB,AMB,AMB,SUB,SUB,SUB,AMB,AMB,AMB,SUB,AMB,
/*c*/  X, SUB,  0,SUB,SUB,AMB,SUB,SUB,AMB,AMB,SUB,AMB,AMB,SUB,AMB,AMB,
/*g*/  X, SUB,SUB,  0,SUB,SUB,AMB,SUB,AMB,SUB,AMB,AMB,SUB,AMB,AMB,AMB,
/*t*/  X, SUB,SUB,SUB,  0,SUB,SUB,AMB,SUB,AMB,AMB,SUB,AMB,AMB,AMB,AMB,
/*m*/  X, AMB,AMB,SUB,SUB,AMB,AMB,AMB,AMB,AMB,SUB,AMB,AMB,AMB,AMB,AMB,
/*r*/  X, AMB,SUB,AMB,SUB,AMB,AMB,AMB,AMB,SUB,AMB,AMB,AMB,AMB,AMB,AMB,
/*w*/  X, AMB,SUB,SUB,AMB,AMB,AMB,AMB,SUB,AMB,AMB,AMB,AMB,AMB,AMB,AMB,
/*s*/  X, SUB,AMB,AMB,SUB,AMB,AMB,SUB,AMB,AMB,AMB,AMB,AMB,AMB,AMB,AMB,
/*y*/  X, SUB,AMB,SUB,AMB,AMB,SUB,AMB,AMB,AMB,AMB,AMB,AMB,AMB,AMB,AMB,
/*k*/  X, SUB,SUB,AMB,AMB,SUB,AMB,AMB,AMB,AMB,AMB,AMB,AMB,AMB,AMB,AMB,
/*v*/  X, AMB,AMB,AMB,SUB,AMB,AMB,AMB,AMB,AMB,AMB,AMB,AMB,AMB,AMB,AMB,
/*h*/  X, AMB,AMB,SUB,AMB,AMB,AMB,AMB,AMB,AMB,AMB,AMB,AMB,AMB,AMB,AMB,
/*d*/  X, AMB,SUB,AMB,AMB,AMB,AMB,AMB,AMB,AMB,AMB,AMB,AMB,AMB,AMB,AMB,
/*b*/  X, SUB,AMB,AMB,AMB,AMB,AMB,AMB,AMB,AMB,AMB,AMB,AMB,AMB,AMB,AMB,
/*n*/  X, AMB,AMB,AMB,AMB,AMB,AMB,AMB,AMB,AMB,AMB,AMB,AMB,AMB,AMB,AMB
};

/* This might be a good idea, non? */
/*    static struct diff_s {
                           int  typ;
                           int  top_pos;
                           int  bot_pos;
                           int  top_char;
                           int  bot_char;
                         } diffs[MAX_SEQ_LENGTH];
*/

                         /*************************/
                         /* Function Declarations */
                         /*************************/


void  gen_alignment( BYTE *backpnt, int m, int n, 
                     char *seq1, char *seq2,
                     int   max_align_size,
                     char *top, char *mid, char *bot, int *l );

void  pq_init( int m, int n );
int   pq_update( int i, int j, int d );
int   pq_remove( int *i, int *j );
void  pq_upheap( int k );
void  pq_downheap( int k );

#define  ind(x,y) ( (x) * (n+1) + (y) )
#define  MAX(x,y) ( ( (x) > (y) ) ? (x) : (y) )

/* Note: seqs and align arrays should be kept below MAX_SEQ_LENGTH  */
/* Hmmm.  maybe better to have stats generated directly from matrix */
/* and then skip align?  Later!                                     */

void lpa_align( char *char_seq1, int m, /* first (top) sequence & length */
                char *char_seq2, int n, /* 2nd (bottom) sequence & length */
                int   max_align_size,   /* max length of top, mid, bot */
                char *top,          
                char *mid, 
                char *bot,
                int  *l                 /* resultant length of top, mid, bot */
              )
                
   {
    int          *d;
    BYTE         *backpnt;
    int           size;
    int           i, j, k;
    int           x, y;
    int           dist;
    int           p;
    int           max_right_lead_offset;
    int           max_down_lead_offset;
    char          *seq1, *seq2;   /* sequential code mapped versions of seqs*/
    int           counter = 0;

    /* allocate storage for matrices */
    size =  (m + 1) * (n + 1);

    d = (int *) calloc_safely( size, sizeof( int ), 
                               "occurred in lpa_align: d" );
    for ( i = 0; i < size; i++ )
        d[i] = -UNSEEN;
    backpnt = (BYTE *) calloc_safely( size, sizeof( BYTE ),
                                      "occurred in lpa_align: backpnt" );

    /* map strings to sequential codes for fast penalty matrix evaluation */
    seq1 = calloc_safely( m, sizeof( char ), "occured in lpa_align: seq1" );
    seq2 = calloc_safely( n, sizeof( char ), "occured in lpa_align: seq2" );
    char_to_sequential( char_seq1, seq1 );
    char_to_sequential( char_seq2, seq2 );

    /* estimate lead offsets */
 
    if ( n > m )
       {
        max_right_lead_offset = n - m;
        max_down_lead_offset = m / 10;
       }
    else
       {
        max_right_lead_offset = n / 10;
        max_down_lead_offset = m - n;
       }

#ifdef DEBUG
    printf( "starting lpa\n" );  fflush( stdout );
#endif
    pq_init( m, n );
    pq_update( 0, 0, UNSEEN );    

    do {
        if ( ! pq_remove( &i, &j ) )
           {
            fprintf( stderr, "error: priority queue empty!\n" );
            exit( 1 );
           }
        if ( (++counter % 10000) == 0 )
            printf( "%d  i,j = (%5d,%5d) dist %8d\n", counter, i, j, dist );
        if ( i < m || j < n )
	   {
            x = ind( i, j );
            d[x] = -d[x];
            if ( d[x] == UNSEEN )
                d[x] = 0;

            if ( j < n )
 	       {
                y = ind( i, j+1 );                        /* RIGHT */
                if ( d[y] < 0 )
		   {
                    if ( i == m )
                        p = 0;
		    else if ( i == 0 && j < max_right_lead_offset )
                        p = 0;
                    else 
                        p = INDEL;
                    dist = d[x] + p;
                    if ( pq_update( i, j+1, dist ) )
	               {
                        d[y] = -dist;
                        backpnt[y] = LEFT;
#ifdef DEBUG
                        printf( "       RIGHT, %d,         %d, %d\n",
                                                   p, dist, d[y] );
#endif
	               }
		   }
                
	       }

            if ( i < m )
	       {
                y = ind( i+1, j );                        /* DOWN */
                if ( d[y] < 0 )
		   {
                    if ( j == n )
                        p = 0;
		    else if ( j == 0 && i < max_down_lead_offset )
                        p = 0;
                    else 
                        p = INDEL;
                    dist = d[x] + p;
                    if ( pq_update( i+1, j, dist ) )
	               {
                        d[y] = -dist;
                        backpnt[y] = UP;
#ifdef DEBUG
                        printf( "       DOWN, %d,         %d, %d\n", 
                                                   p, dist, d[y] );
#endif
	               }
		   }
	       }

            if ( i < m && j < n )
	       {
                y = ind( i+1, j+1 );                      /* DIAG */
                if ( d[y] < 0 )
		   {
     	            /*
                    if ( seq1[i] == seq2[j] )
                        p = 0;
                    else if ( ISAMB( seq1[i] ) || ISAMB( seq2[j] ) )
                        p = AMB_SUBST;
                    else
                        p = SUBST_ERR;
    		    */
                    p = penalty[seq1[i]][seq2[j]];
                    dist = d[x] + p;
                    if ( pq_update( i+1, j+1, dist ) )
	               {
                        d[y] = -dist;
                        backpnt[y] = DIAG;
#ifdef DEBUG
                        printf( "       DIAG, %d,           %d, %d\n", 
                                p, dist, d[y] );
#endif
	               }
		   }
	       }
	   }
       }
    while ( i < m || j < n );

#ifdef DEBUG
    printf( "Done lpa!\n" );
#endif

    gen_alignment( backpnt, m, n, char_seq1, char_seq2, max_align_size, 
                   top, mid, bot, l );
#ifdef DEBUG
    printf( "size is %d, align l is %d\n", size, *l );
#endif

    free( d );
    free( backpnt );
    free( seq1 );
    free( seq2 );
   }



void  gen_alignment( BYTE *backpnt, int m, int n, 
                     char *seq1, char *seq2,
                     int   max_align_size,
                     char *top, char *mid, char *bot, int *l )
   {
    int  i, j, k;
    char *s1, *s2;
    
    i = m;
    j = n;
    k = 0;
    s1 = seq1 + m;
    s2 = seq2 + n;
    while ( k < max_align_size && (i > 0 || j > 0) )
       {
#ifdef DEBUG
        printf( "gen al: %d %d: %d\n", i, j, backpnt[ ind(i,j) ] );
#endif
        switch ( backpnt[ ind(i,j) ] )
           {
            case  UP:
                          top[k] = *--s1;
                          mid[k] = INDEL_IND;
                          bot[k] = INDEL_SPACER;
                          i--;
                          break;
            case  LEFT:
                          top[k] = INDEL_SPACER;
                          mid[k] = INDEL_IND;
                          bot[k] = *--s2;
                          j--;
                          break;
            case  DIAG:
                          top[k] = *--s1;
                          bot[k] = *--s2;
                          if ( *s1 == *s2 )
                              mid[k] = OK_IND;
                          else if ( ISAMB(SEQUENT(*s1)) || ISAMB(SEQUENT(*s2)))
                              mid[k] = AMB_IND;
                          else
                              mid[k] = WRONG_IND;
                          i--; j--;
                          break;
	   }
        k++;
       }
    if ( k >= max_align_size )
       {
        fprintf( stderr, "%s (gen_alignment): alignment array overflow!\n",
                         progname );
        exit( 1 );
       }
    *l = k;
   }

#define SUBST_TYPE      1
#define TOP_INS_TYPE    2
#define BOT_INS_TYPE    3

void  get_align_stats( char *top, char *mid, char *bot, int l, 
                       int m, int n,
                       ALIGN_STATS *s )

   {
    int  i, k;
    int  a, b;
    int  top_pos, bot_pos;

    s->top_end_off = s->bot_end_off = s->top_beg_off = s->bot_beg_off = 0;
    s->num_ok = s->num_wrong = s->num_indels = s->num_ambs = 0;
    s->num_seq1_inserts = s->num_seq2_inserts = 0;
    s->num_seq1_ambs = s->num_seq2_ambs = 0;

    /* advance past beginning offsets, and count them.  Don't forget that */
    /* top, mid & bot are backwards!                                      */
    l--;  
    while ( l >= 0 && mid[l] == INDEL_IND )
        if ( top[l--] == INDEL_SPACER )
            s->top_beg_off++;
        else
            s->bot_beg_off++;
    b = l;    

    /* l now points at first "interior" char */
    
    i = 0;
    while ( i <= l && mid[i] == INDEL_IND )
        if ( top[i++] == INDEL_SPACER )
            s->top_end_off++;
        else
            s->bot_end_off++;
    a = i;

    /* i now points at last "interior" character */

    top_pos = s->bot_beg_off + 1;  /* note the "switched" sense here */
    bot_pos = s->top_beg_off + 1;

    while ( l >= i )
       {
        if ( mid[l] != OK_IND )
           {
#ifdef DEBUG         /* This should be put into a separate routine! */
            printf( "%c: %d %c - %d %c\n", mid[l], top_pos, top[l], 
                                                   bot_pos, bot[l] );
#endif
            if ( top[l] != INDEL_SPACER )
                top_pos++;
            if ( bot[l] != INDEL_SPACER )
                bot_pos++;
           }
        else
           {
            top_pos++;
            bot_pos++;
           }
        l--;
       }
#ifdef DEBUG
    printf( "i = %d, l = %d\n", i, l );
    printf( "a = %d, b = %d\n", a, b );
    printf( "%c %c %c   %c %c %c\n", top[a], mid[a], bot[a], 
                                     top[b], mid[b], bot[b] );
#endif
    for ( k = a; k < b; k++ )
        switch( mid[k] )
	   {
	    case  OK_IND:      s->num_ok++;
                               break;
            case  INDEL_IND:   s->num_indels++;
                               if ( top[k] == INDEL_SPACER )
                                   s->num_seq1_inserts++;
                               else
                                   s->num_seq2_inserts++;
                               break;
            case  AMB_IND:     s->num_ambs++;
                               if ( ISAMB( top[k] ) )
                                   s->num_seq1_ambs++;
                               else
                                   s->num_seq2_ambs++;
                               break;
            case  WRONG_IND:   s->num_wrong++;
                               break;

	   }
   }


void  print_align_stats( ALIGN_STATS *a )
   {
    printf( "seq1 offsets:            %4d    %4d\n", a->top_beg_off, 
                                                     a->top_end_off );
    printf( "seq2 offsets:            %4d    %4d\n", a->bot_beg_off, 
                                                     a->bot_end_off );
    printf( "Matches:                 %4d\n", a->num_ok );
    printf( "Mismatches:              %4d\n", a->num_indels + a->num_ambs
                                              + a->num_wrong );
    printf( "Errors:                  %4d\n", a->num_indels + a->num_wrong );
    printf( "Substitution Errors:     %4d\n", a->num_wrong );
    printf( "Indels:                  %4d\n", a->num_indels );
    printf( "Insertions in seq 1:     %4d\n", a->num_seq1_inserts );
    printf( "Insertions in seq 2:     %4d\n", a->num_seq2_inserts );
    printf( "Ambiguity substitutions: %4d\n", a->num_ambs );
    printf( "Ambiguities in seq 1:    %4d\n", a->num_seq1_ambs );
    printf( "Ambiguities in seq 2:    %4d\n", a->num_seq2_ambs );
   }


/* Prints out the alignment that was constructed backwards in top, */
/* mid & bot starting at position n-1                              */

char empty[] = { 0 };   /* this is acting up - better be careful! */

void  print_alignment( FILE *f, char *top_name, char *bot_name,
                       char *top, char *mid, char *bot, int n )

   {
    int  i, j, k, l;

    i = j = k = n;
    do {
        fprintf( f, "%-18.18s ", top_name );
        l = 0;
        while ( i > 0 && l++ < 60 )
            putc( top[--i], f );
        putc( '\n', f );

        fprintf( f, "%-18.18s ", empty );
        l = 0;
        while ( j > 0 && l++ < 60 )
            putc( mid[--j], f );
        putc( '\n', f );

        fprintf( f, "%-18.18s ", bot_name );
        l = 0;
        while ( k > 0 && l++ < 60 )
            putc( bot[--k], f );
        putc( '\n', f );
        for ( l = 0; l < 79; l++ )
            putc( '_', f );
        putc( '\n', f );
       }
    while ( i > 0 );
   }

/* convention: bottom is raw, top is corrected sequence.                     */
/*    meanings:  ins - ABI made an insertion _error_                         */
/*               del - ABI made a deletion _error_                           */
/*               sub - ABI made a substition error (even if top is N)        */

void  print_diffs( ALIGN_STATS *a, char *top, char *mid, char *bot, int n )

   {
    int  i;       /* position counter for top (corrected) */
    int  j;       /* position counter for bot (corrected) */
    int  n_stop;  /* end of overlap region                */

    /* remember, now, that top, mid & bot are backwards! */

    i = j = 1;
    if ( a->top_beg_off > 0 )
       {
        n -= a->top_beg_off;
        j += a->top_beg_off;
       }
    else if ( a->bot_beg_off > 0 )
       {
        n -= a->bot_beg_off;
        i += a->bot_beg_off;
       }
    n_stop = MAX( a->top_end_off, a->bot_end_off );

    while ( --n >= n_stop )
       {
        switch( mid[n] )
	   {
	    case  OK_IND:      i++; j++;
                               break;
            case  INDEL_IND:   if ( top[n] == INDEL_SPACER )
                                   printf( "%4d %4d ins - %c\n", 
                                            i, j++, bot[n] );
                               else
                                   printf( "%4d %4d del %c -\n", 
                                            i++, j, top[n]  );
                               break;
            case  AMB_IND:     printf( "%4d %4d amb %c %c\n", 
                                       i++, j++, top[n], bot[n] );
                               break;
            case  WRONG_IND:   printf( "%4d %4d sub %c %c\n", 
                                       i++, j++, top[n], bot[n] );
                               break;

	   }
       }
   }



                        /***************************/
                        /* Priority Queue Routines */
                        /***************************/


#define PQ_TREE_SIZE   64000
#define MAX_PRIOR         -1

int     pq_n = 0;
int     pq_tree_dist[PQ_TREE_SIZE+1];
int     pq_tree_x[PQ_TREE_SIZE+1];
int    *pq_lkup;
short  *pq_i;
short  *pq_j;
int     pq_mat_n;
int     pq_mat_m;
int     pq_mat_size;



#define pq_ind(x,y) ((x)*(pq_mat_n+1) + (y))

void pq_init( int m, int n )

   {
    pq_mat_m = m;
    pq_mat_n = n;
    pq_mat_size = (n + 1) * (m + 1);
#ifdef DEBUG
    printf( "pq_mat size is %d\n", pq_mat_size );
#endif
    pq_lkup = (int *) calloc_safely( pq_mat_size, sizeof(int),
                                     "occurred in pq_ind: pq_lkup" );
    pq_i = (short *) calloc_safely( pq_mat_size, sizeof(short),
                                     "occurred in pq_ind: pq_i" );
    pq_j = (short *) calloc_safely( pq_mat_size, sizeof(short),
                                     "occurred in pq_ind: pq_j" );
    pq_n = 0;
   }

/* returns 1 if an insert is made or if the priority is raised     */
/* (if dist is < the current value), and returns 0 if no change is */
/* made.                                                           */

int pq_update( int i, int j, int dist )

   {
    int x, k;

    x = pq_ind(i,j);
    /* Is  i,j (x rather) in pq already? */
    if ( (k = pq_lkup[x]) > 0 )
       {
        if ( pq_tree_dist[k] <= dist )
            return( 0 );
        pq_tree_dist[k] = dist;
        pq_upheap( k );
       }
    else
       {
        if ( pq_n >= PQ_TREE_SIZE )
           {
            fprintf( stderr, "pq tree overflow, i = %d, j = %d, dist = %d\n",
                              i, j, dist );
            exit( 1 );
           }
        x = pq_ind(i,j);
        pq_tree_dist[++pq_n] = dist;
        pq_tree_x[pq_n] = x;
        pq_lkup[x] = pq_n;
        pq_i[x] = i;
        pq_j[x] = j;
        pq_upheap( pq_n );
       }

    return( 1 );
   }


int pq_remove( int *i, int *j )

   {
    int x;
    
    if ( pq_n < 1 )
        return( 0 );
    x = pq_tree_x[1];
    *i = pq_i[x];
    *j = pq_j[x];
    pq_lkup[x] = 0;   /* indicates that i,j is no longer in pq */
                      /* - used by pq_update().                */
    pq_tree_dist[1] = pq_tree_dist[pq_n];  /* now put last guy at top */
    pq_tree_x[1]    = pq_tree_x[pq_n];     
    pq_lkup[pq_tree_x[pq_n]] = 1;
    pq_n--;
    pq_downheap( 1 );                      /* and bubble him down */
    return( 1 );
   }

   
void  pq_upheap( int k )
   {
    int d, j, x;
    
    d = pq_tree_dist[k];
    x = pq_tree_x[k];
    pq_tree_dist[0] = MAX_PRIOR; 
    while ( pq_tree_dist[ j = ( k / 2 ) ] > d )
       {
        pq_lkup[ pq_tree_x[j] ] = k;
        pq_tree_dist[k] = pq_tree_dist[j];
        pq_tree_x[k] = pq_tree_x[j];
        k = j;
       }
    pq_lkup[x] = k;
    pq_tree_dist[k] = d;
    pq_tree_x[k] = x;
   }

void  pq_downheap( int k )
   {
    int d, x, j;

    d = pq_tree_dist[k];
    x = pq_tree_x[k];
    while ( k <= pq_n / 2 )
       {
        j = k + k;            /* k <<1 */
        if ( j < pq_n )
            if ( pq_tree_dist[j] > pq_tree_dist[j+1] )
                j++;
        if ( d < pq_tree_dist[j] ) 
            break;
        pq_lkup[ pq_tree_x[j] ] = k;
        pq_tree_dist[k] = pq_tree_dist[j];
        pq_tree_x[k] = pq_tree_x[j];
        k = j;
       }
     pq_tree_dist[k] = d;
     pq_tree_x[k] = x;
     pq_lkup[x] = k;
   }

#ifdef DEBUG
void pq_dump()

   {
    int i, j;
    
    printf( "pq tree size: %d\n", pq_n );
    for ( i = 0; i <= pq_n; i++ )
        printf( "%3d: %4d  (%3d %3d)\n", 
                 i, pq_tree_dist[i], pq_i[pq_tree_x[i]], 
                 pq_j[pq_tree_x[i]]);
#ifdef SUPERDEBUG
    for ( i = 1; i <= pq_mat_m; i++ )
        for ( j = 1; j <= pq_mat_n; j++ )
            printf( "(%2d %2d) [%2d] %2d (%2d %2d)\n", 
                     i, j, pq_ind(i,j), pq_lkup[pq_ind(i,j)],
                     pq_i[pq_ind(i,j)], pq_j[pq_ind(i,j)] );
#endif
   }
#endif



