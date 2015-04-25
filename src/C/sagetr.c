/* Program:      sagetr.c                                                   */
/* Programmer:   Sean R. McCorkle                                           */
/* Language:     C                                                          */
/*                                                                          */
/* Description:  Provide statistics on SAGE tages occuring in specified     */
/*               gene (mRNA or cDNA data)                                   */
/*                                                                          */
/* Usage:                                                                   */
/*               sagetr [options] [<mRNA file> ...]                         */
/*                                                                          */
/*               where <mRNA file> is in FASTA format.  Stdin is scanned    */
/*               if - is specified or no files are specfied                 */
/*                                                                          */
/*               Options:                                                   */
/*                                                                          */
/*                 -e <n>  tag extent (after enzyme sequence, default 10)   */
/*                 -p <seq> puncutating enzyme sequence (deault "CATG")     */
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
/*                 -T       dump out trie - copious output for debuging     */
/*                                                                          */
/* $Id: sagetr.c,v 0.3 2001/04/10 17:17:08 mccorkle Exp mccorkle $          */
/****************************************************************************/


static char sagetr_rcs_id[] =
      "$Id: sagetr.c,v 0.3 2001/04/10 17:17:08 mccorkle Exp mccorkle $";

#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifndef SOLARIS        /* Sun's are really a pain in the neck sometimes */
#include <unistd.h>
#endif
#include "seqlib.h"
#include "trie.h"


#define  MAX_SEQ_LEN   200000
#define  MAX_SEQUENCES 100000

typedef enum { FORWARD = 'F', REVERSE = 'R' } DIRECTION;
typedef enum { TRUE = 1, FALSE = 0 } BOOLEAN;
 
typedef struct {
                 int             string_id;
                 int             pos;
                 DIRECTION       dir;
                 BOOLEAN         real_tag;  /* F => not the real tag */
               } POS_REC;

/* globals which are potentially set by command line options */

char  *re_seq         = "CATG";  /* -p */
int    tag_extent     = 10;
int    list_hits      = 0;       /* -c */
int    list_negs      = 0;       /* -n */
int    verbose        = 0;       /* -v */
int    print_suffixes = 0;       /* -s */
int    print_tags     = 0;       /* -t */
int    print_pos      = 0;       /* -l */
int    summarize      = 0;       /* -S */
int    dump_trie      = 0;       /* -T */
int    n_files;
char **files;

/* globals for calculations */

int    n_tags = 0;
int    n_mult_tags = 0;
int    n_uniq_tags = 0;
int    n_really_uniq_tags = 0;


char *seqnames[MAX_SEQUENCES];

POS_REC  *pos_rec( int id, int pos, DIRECTION d, BOOLEAN real_tag )
   {
    POS_REC *p;

    p = (POS_REC *) malloc_safely( sizeof( POS_REC ), "pos_rec" );
    p->string_id = id;
    p->pos = pos;
    p->dir = d;
    p->real_tag = real_tag;
    return( p );
   }


void  print_pos( void *d )
   {
    POS_REC *p = (POS_REC *) d;
    printf( "  (%d %s %d %c %d)\n", p->string_id, seqnames[p->string_id],
                                     p->pos, p->dir, p->real_tag );
   }


void  analyze( TR_LEAF *l )
   {
    int      n;
    int      n_real;
    TR_LEAF *q;

    if ( l == NULL )
        return;
    for (  q = l, n = 0, n_real = 0 ;  q != NULL;  q = q->next, n++  )
        if ( ((POS_REC *)(q->data))->real_tag )
            n_real++;

    if ( n_real > 0 )
        n_tags += n_real;

    if ( n_real == 1 )
       {
        n_uniq_tags++;
        if ( n == 1 )
            n_really_uniq_tags++;
       }
    else if ( n_real > 1 )
       {
        n_mult_tags++;
       }
    if ( n_real >= 1 && print_tags )
       {
        printf( "%-14.14s %10d %10d\n", DFSpath(), n_real, n - n_real );
        if ( verbose )
           {
            for ( q = l; q != NULL; q = q->next )
               print_pos( q->data );
            printf( "----\n" );
           }
       }
   }


void  parse_args( int argc, char **argv )
   {
    int          c;
    int          i;
    extern char *optarg;
    extern int   optind;
          
    while ( (c = getopt(argc, argv, "ce:lnp:hstSTv") ) != -1 )
        switch (c) 
           {
            case 'c':
                       list_hits = 1;
                       break;
            case 'e':
                       tag_extent = atoi( optarg );
                       break;
            case 't': 
                       print_pos = 1;
                       break;
            case 'n':
                       list_negs = 1;
                       break;
            case 'p':
                       re_seq = strdup( optarg );
                       uppercase( re_seq, re_seq );
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
            case 'T': 
                       dump_trie = 1;
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


int  get_seq( int seqf, char *seq, int *len, char *hdr )
   {
    int   ret;
    char *b;

    ret = read_seq( seqf, MAX_HDR_LENGTH, hdr, MAX_SEQ_LEN, seq );
    if ( ret > 0 ) 
       {
        *len = strlen( seq );
        uppercase( seq, seq );
#ifdef DEBUG
        if ( b = index( hdr, ' ' ) )
  	    *b = '\0';                 /* cut off name at first blank */
#endif
       }
    return( ret );
   }


int  main( int argc, char **argv )
   {
    TR_NODE       *t;
    int            seqf;
    int            i;
    static char    seq[MAX_SEQ_LEN];
    int            len;    /* length of sequence being processed */
    int            re_len; /* length of punctuating enzyme sequence */
    int            id;     /* sequence id (incrementing counter) */
    char          *p;      /* position of punc. restriction site in seq */
    static char    tag[200];
    int            tag_len;
    int            negative_fs = 0;  /* count of seqs with no re site forward*/
    POS_REC       *r;
    static char    seqnam[MAX_HDR_LENGTH];
    static char   *poly_a = "AAAAAAAAAAAAAAAAAAAA";

    parse_args( argc, argv );
    re_len = strlen( re_seq );
    tag_len = re_len + tag_extent;
    t = create_trie();
    id = 0;

    /* read all sequences, scan for enzyme sites, add to trie where found */

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
                if ( list_hits )
                    printf( ">Y %s\n", seqnam );
                do {
                    strncpy( tag, p + re_len, tag_extent );
#ifdef EXTREME_DEBUG
                    printf( "id %d, tag is [%s]\n", id, tag);
#endif
                    r = pos_rec( id, (int)(p-seq) + 1, FORWARD, FALSE );
                    add_trie_string( t, tag, r );
                   }
                while ( p = strstr( p + 1, re_seq ) );
                r->real_tag = TRUE;
                if ( print_suffixes )
                   {
                    printf( ">tag-tail: %s\n", seqnames[id] );
                    print_seq( stdout, seq, r->pos, len, 0, 0, 0, 0 );
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

    /* All sequences now read, and the trie is built */

    if ( dump_trie )
        print_trie( t, 0, print_pos );  

    if ( print_tags )
       {
        printf( "\n\n%s+%d tags\n\n", re_seq, tag_extent );
        printf( "                   occurs     occurs\n" );
        printf( "tag                as tag     nontag\n" );
        printf( "--------------     ------     ------\n" );
       }

    DFSprocess( t, analyze );

    if ( print_tags )
        printf( "\n" );

    if ( summarize )
       {
        printf( "%s+%d: ", re_seq, tag_extent );
        printf( 
          "%d sequences, %d tags, %d negatives,  %d unique, %d multiples\n",
           id,            n_tags, negative_fs, n_uniq_tags, n_mult_tags
              );  
        printf( "        (%d highly unique tags)\n", n_really_uniq_tags );
       }
   }





