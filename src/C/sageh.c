/* Program:      sageh.c                                                    */
/* Programmer:   Sean R. McCorkle                                           */
/* Language:     C                                                          */
/*                                                                          */
/* Description:  Provide statistics on SAGE tages occuring in specified     */
/*               gene (mRNA or cDNA data)                                   */
/*                                                                          */
/* Usage:                                                                   */
/*               sageh [options] [<mRNA file> ...]                          */
/*                                                                          */
/*               where <mRNA file> is in FASTA format.  Stdin is scanned    */
/*               if - is specified or no files are specfied                 */
/*                                                                          */
/*               Options:                                                   */
/*                                                                          */
/*                 -a <seq> anchoring enzyme sequence (deault "CATG")       */
/*                 -e <n>  tag extent (after enzyme sequence, default 10)   */
/*                                                                          */
/*               Output choices:                                            */
/*                                                                          */
/*                 -c       list seqnames which do contain ensymze site     */
/*                 -l       print tags and locations                        */
/*                 -n       list seqnames which do NOT contain enzyme site  */
/*                 -s       print_suffixes = 1;                             */
/*                 -t       print_tags and occurence rates                  */
/*                 -v       verbose mode                                    */
/*                 -S       print summary (stats)                           */
/*                 -D       dump out hash table (debugging)                 */
/*                                                                          */
/*                                                                          */
/* $Id: sageh.c,v 1.3 2001/05/23 18:03:05 mccorkle Exp mccorkle $ */
/****************************************************************************/

#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifndef SOLARIS        /* Sun's are really a pain in the neck sometimes */
#include <unistd.h>
#endif
#ifdef  IRIX6
#include <strings.h>
#endif
#include <search.h>
#include "seqlib.h"

static char sageh_rcs_id[] =
      "$Id: sageh.c,v 1.3 2001/05/23 18:03:05 mccorkle Exp mccorkle $";

#define  MAX_SEQ_LEN   200000
#define  MAX_SEQUENCES 100000
#define  MAX_TAGS      (5 * MAX_SEQUENCES)


typedef enum { FORWARD = 'F', REVERSE = 'R' } DIRECTION;
typedef enum { TRUE = 1, FALSE = 0 } BOOLEAN;

/* POS_REC - position record: one for each tag "hit" anywhere in the data*/

typedef struct p_rec {
                      int             string_id;
                      int             pos;
                      DIRECTION       dir;
                      BOOLEAN         real_tag;  /* F => not the real tag */
                      struct p_rec   *next;
                     } POS_REC;

/* TAG_REC - one for each tag sequence in hash table.  This contains a link */
/* to a list of POS_RECs (above).  The routine count_tags() fills the       */
/* counters in each POS_REC                                                 */

typedef struct t_rec {
                      int             n_occurs;  /* times it occurs in seqs  */
                      int             n_real;    /* times its a real tag     */
                      POS_REC        *plist;     /* links to list of pos recs*/
                     } TAG_REC;


/* globals which are potentially set by command line options */

char  *re_seq          = "CATG";  /* -a */
int    tag_extent      = 10;
int    list_hits       = 0;       /* -c */
int    list_negs       = 0;       /* -n */
int    verbose         = 0;       /* -v */
int    print_suffixes  = 0;       /* -s */
int    print_tags      = 0;       /* -t */
int    print_positions = 0;       /* -l */
int    summarize       = 0;       /* -S */
int    dump_hash_table = 0;       /* -D */
int    short_names     = 1;       /* ? */
char   name_seperator  = ' ';     /* ? */
int    n_files;
char **files;

/* globals for calculations */

int    n_tags             = 0;
int    n_mult_tags        = 0;
int    n_uniq_tags        = 0;
int    n_really_uniq_tags = 0;
int    n_positives        = 0;
int    negative_fs        = 0;    /* count of seqs with no re site forward*/

int    n_sequences        = 0;
char  *seqnames[MAX_SEQUENCES];

int    curr_tag           = 0;   /* redundant with n_tags. Fix! */
char  *tags[MAX_TAGS];

void  bailout( char *msg )
   {
    fprintf( stderr, "%s\n", msg ); 
    exit( 1 );
   }

POS_REC  *pos_rec( int id, int pos, DIRECTION d, BOOLEAN real_tag )
   {
    POS_REC *p;

    p = (POS_REC *) malloc_safely( sizeof( POS_REC ), "pos_rec" );
    p->string_id = id;
    p->pos = pos;
    p->dir = d;
    p->real_tag = real_tag;
    p->next = NULL;
    return( p );
   }

TAG_REC  *tag_rec( POS_REC *r )
   {
    TAG_REC *t;

    t = (TAG_REC *) malloc_safely( sizeof( TAG_REC ), "tag_rec" );
    t->n_occurs = 0;
    t->n_real = 0;
    t->plist = r;
    return( t );
   }


/*void  print_pos( void *d )*/

void  print_pos( POS_REC *p )
   {
     /*    POS_REC *p = (POS_REC *) d; */
    printf( "  (%d %s %d %c %d)\n", p->string_id, seqnames[p->string_id],
                                     p->pos, p->dir, p->real_tag );
   }

void print_poslist( POS_REC *r )
   {
    for ( ; r != NULL; r = r->next )
        print_pos( r );
   }

/* for sorting alphabetically */

int  alpha_cmp( const void *a, const void *b )
   {
    return( strcmp( *((char **) a), *((char **) b) ) );
   }

int  get_seq( int seqf, char *seq, int *len, char *hdr )
   {
    int   ret;
    char *b;

    ret = read_seq( seqf, MAX_HDR_LENGTH, hdr, MAX_SEQ_LEN, seq );
    if ( ret > 0 ) 
       {
        *len = strlen( seq );
        uppercase( seq, seq );
        if ( short_names )
            if ( b = index( hdr, ' ' ) )
                *b = '\0';                 /* cut off name at first blank */
       }
    return( ret );
   }


void  parse_args( int argc, char **argv )
   {
    int          c;
    int          i;
    extern char *optarg;
    extern int   optind;
          
    while ( (c = getopt(argc, argv, "a:cDe:lnhstSv") ) != -1 )
        switch (c) 
           {
            case 'a':
                       re_seq = strdup( optarg );
                       uppercase( re_seq, re_seq );
                       break;
            case 'c':
                       list_hits = 1;
                       break;
            case 'D': 
                       dump_hash_table = 1;
                       break;
            case 'e':
                       tag_extent = atoi( optarg );
                       break;
            case 'l': 
                       print_positions = 1;
                       break;
            case 'n':
                       list_negs = 1;
                       break;
            case 'h': 
                       fprintf( stderr, "no help yet. sorry.\n" );
                       exit( 0 ); 
                       break;
            case 's': 
                       print_suffixes = 1;
                       break;
            case 'S': 
                       summarize = 1;
                       break;
            case 't': 
                       print_tags = 1;
                       break;
            case 'v': 
                       verbose = 1;
                       break;
            case '?': 
                       fprintf( stderr, "Goodbye, Mr. Anderson\n" );
                       exit( 0 );
           }
    n_files = argc - optind;    
    files = (char **) malloc( n_files * sizeof( char * ) );
    i = 0;
    while( optind < argc )
        files[i++] = strdup( argv[optind++] ); 
   }


void  initialize_hash_table( void )
   {
    if ( ! hcreate( MAX_TAGS ) )
       {
        perror( "could not create hash table\n" );
        exit( errno );
       }
   }

/* Lookup tag string in hash table, return pointer to TAG_REC if found,    */
/* otherwise return null                                                   */

TAG_REC *lookup_tag( char *tag )
   {
    ENTRY    ent, *e;

    ent.key = tag;                           /* load entry for initial srch*/
    ent.data = NULL;
    if ( e = hsearch( ent, FIND ) )
        return( (TAG_REC *)(e->data) );
    else
        return( (TAG_REC *) NULL );
   }


void  enter_tag( char *tag, POS_REC *r )
   {
    ENTRY    ent, *e;
    TAG_REC  *t;

#ifdef EXTREME_DEBUG
    printf( "enter tag [%s]\n", tag );
    print_pos( r );
#endif

    if ( t = lookup_tag( tag ) )           /* is the tag in the hash?*/
       {                           
#ifdef EXTREME_DEBUG
        printf( "#############Tag Found!!!!!\n" );
#endif
        r->next = t->plist;                   /* if so, add new pos_rec to */
        t->plist = r;                         /* the existing entry list */
#ifdef EXTREME_DEBUG
        print_poslist( t->plist );
#endif
       }
    else 
       {
#ifdef EXTREME_DEBUG
        printf( "Not found.  Inserting\n" );
#endif
        if ( curr_tag >= MAX_TAGS )
           {
            fprintf( stderr, "number of unique tags exceeded %d\n", MAX_TAGS );
            exit( 1 );
           }
        ent.key = tags[curr_tag++] = strdup( tag );
        ent.data = (char *) tag_rec( r );
        if ( ! hsearch( ent, ENTER ) )
           {
            perror( "Can't enter tag, hsearch( tag, ENTER)" );
            exit( errno );
           } 
       }
   }

void  dump_data( void )
   {
    int      i;
    TAG_REC *t;

    printf( "Curr tag %d\n", curr_tag );
    for ( i = 0; i < curr_tag; i++  )
       {
        printf( "tag: [%s] %d\n", tags[i], i );
        if ( !(t = lookup_tag( tags[i] ) ) )
            bailout( "dump_data: couldn't retrieve key" );
        print_poslist( t->plist );
       }
   }


void terminate_hash_table( void )
   {
    hdestroy();
   }

void  count_tags( void )
   {
    TAG_REC *t;
    POS_REC *p;
    int      i;

    if ( verbose )
        printf( "Counting...\n" );  fflush( stdout );

    for ( i = 0; i < curr_tag; i++  )                /* for each tag */
       {
        if ( !(t = lookup_tag( tags[i] ) ) )
            bailout( "count_tags: could not retrieve key" );
        t->n_occurs = t->n_real = 0;
        for ( p = t->plist; p != NULL; p = p->next ) /* for each position rec*/
           {
            t->n_occurs++;
            if ( p->real_tag )
                t->n_real++;
           }

        /* Update global counters */

        if ( t->n_real > 0 )
            n_tags++;
        if ( t->n_real == 1 )
           {
            n_uniq_tags++;
            if ( t->n_occurs == 1 )
                n_really_uniq_tags++;
           }
        else if (  t->n_real > 1 )
            n_mult_tags++;
       }
    if ( verbose )
        printf( "...Finished!\n" );  fflush( stdout );
   }


void  sort_tags( void )
   {
    if ( verbose )
        printf( "Starting qsort...\n" );
    qsort( tags, curr_tag, sizeof(char *), alpha_cmp );
    if ( verbose )
        printf( "...Finished qsort\n" );
   }



void  print_tag_report( void )
   {
    TAG_REC *t;
    int      i;

    printf( "\n\n%s+%d tags\n\n", re_seq, tag_extent );
    printf( "                   occurs     occurs\n" );
    printf( "tag                as tag     nontag\n" );
    printf( "--------------     ------     ------\n" );

    for ( i = 0; i < curr_tag; i++  )           /* for each tag */
       {
        if ( !(t = lookup_tag( tags[i] ) ) )
            bailout( "print_tag_report: could not retrieve key" );
        if ( t->n_real > 0 )
            printf( "%-14.14s %10d %10d\n", tags[i], t->n_real, 
                                            t->n_occurs - t->n_real );
       }
    printf( "\n" );
   }

void  print_position_report( void )
   {
    TAG_REC *t;
    POS_REC *p;
    int      i;

    printf( "seq               tag                pos   \n" );
    printf( "--------------    --------------     ------\n" );

    for ( i = 0; i < curr_tag; i++  )           /* for each tag */
       {
        if ( !(t = lookup_tag( tags[i] ) ) )
            bailout( "print_tag_report: could not retrieve key" );
        if ( t->n_real > 0 )
           {
            for ( p = t->plist; p != NULL; p = p->next )
               {
                if ( p->real_tag )
                    printf( "%-14.14s    %-14.14s %10d\n", 
                             seqnames[p->string_id], tags[i], p->pos );
               }
           }
       }
    printf( "\n" );
   }


void  print_summary( void )
   {
    printf( "%s+%d: ", re_seq, tag_extent );
    printf( 
          "%d sequences, %d tags, %d negatives,  %d unique, %d multiples\n",
           n_sequences,   n_positives, negative_fs, n_uniq_tags, n_mult_tags
              );  
    printf( "        (%d highly unique tags)\n", n_really_uniq_tags );
   }


int  main( int argc, char **argv )
   {
    int            seqf;
    int            i;
    static char    seq[MAX_SEQ_LEN];
    int            len;    /* length of sequence being processed */
    int            re_len; /* length of punctuating enzyme sequence */
    int            id;     /* sequence id (incrementing counter) */
    char          *p;      /* position of punc. restriction site in seq */
    static char    tag[200];
    int            tag_len;
    static char    seqnam[MAX_HDR_LENGTH];
    static char   *poly_a = "AAAAAAAAAAAAAAAAAAAA";
    POS_REC       *r;

    parse_args( argc, argv );
    re_len = strlen( re_seq );
    tag_len = re_len + tag_extent;

    initialize_hash_table();

    id = 0;

    /* read all sequences, scan for enzyme sites, add to trie where found */

    if ( verbose )
        printf( "reading sequences...\n" );
    for ( i = 0; i < n_files; i++ )    /* for each file */
       {
        if ( ( seqf = open_seqf( files[i], STD_ALLOWED, STD_IGNORED ) < 0 ) )
           {
            fprintf( stderr, "%s: fatal; ", progname );
            perror( files[i] );
            exit( errno );
           }

        while ( get_seq( seqf, seq, &len, seqnam ) > 0 )   /* for each seq */
           {
            if ( id >= MAX_SEQUENCES )
               {
                fprintf( stderr, "%s: fatal; excceded max of %d sequences\n",
                                   progname, MAX_SEQUENCES );
                exit( errno );
               }
            strcat( seq, poly_a ); /* fix this! */
#ifdef EXTREME_DEBUG
            printf( "sequence loaded\n" );
            print_seq( stdout, seq, 0, 0, 5, 0, 0, 0 );
#endif
            seqnames[++id] = strdup( seqnam );
            r = NULL;
            if ( p = strstr( seq, re_seq ) )  
	       {                                /* enzyme site found */
                n_positives++;
                if ( list_hits )
                    printf( ">Y %s\n", seqnam );
         
                /* for each enzyme site (not just the 3'-most) */

                do {
                    strncpy( tag, p + re_len, tag_extent );
#ifdef EXTREME_DEBUG
                    printf( "id %d, tag is [%s]\n", id, tag);
#endif
                    r = pos_rec( id, (int)(p-seq) + 1, FORWARD, FALSE );

                    enter_tag( tag, r );    /* enter into hash table */
                   }
                while ( p = strstr( p + 1, re_seq ) );

                r->real_tag = TRUE;    /* adjustment: last one WAS a real tag*/
                if ( print_suffixes )
                   {
                    printf( ">tag-tail: %s\n", seqnames[id] );
                    print_seq( stdout, seq, r->pos + re_len - 1, len, 
                                            0, 0, 0, 0 );
                   }
               }
            else                                /* no enzyme site found! */
               {
                negative_fs++;
                if ( list_negs )
                    printf( ">N %s\n", seqnam );
               }
           }
        close_seqf( seqf );
       }
    if ( verbose )
        printf( "...finished reading sequences\n" );

    n_sequences = id;

    if ( dump_hash_table )
        dump_data();

    count_tags();

    sort_tags();

    if ( print_tags )
        print_tag_report();

    if ( print_positions )
        print_position_report();

    if ( summarize )
        print_summary();

    terminate_hash_table();

    if ( verbose )
        printf( "Program ends\n" );
   }
