/* Program:     restr                                                        */
/* Programmer:  Sean R. McCorkle                                             */
/* Language:    C                                                            */
/* Description: find restriction enzyme sites in DNA sequence                */
/* Usage:       restr [-fsr] <seq file>                                      */
/*                                                                           */
/* TODO:  comments in restr_enzyme file                                      */
/*        check: re length > MAX_RESTR_LEN?                                  */
/*        restr_enzyme file in /usr/local/seq/lib                            */
/*        option for override restr_enzyme file                              */
/*        N match in restr enzyme file                                       */
/*        option for center position                                         */
/*                                                                           */
/*                                                                           */
/* $Id$ */
/*****************************************************************************/

#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifndef SOLARIS        /* Sun's are really a pain in the neck sometimes */
#include <unistd.h>
#endif
#include "seqlib.h"
#include "trie.h"

#define  MAX_SEQ_LEN               5000000    /* max input sequence length */
#define  MAX_FILE_LEN                  256    /* max length of file name */
/*#define  MAX_RESTR_LEN                  10 */   /* max length of restr enzyme */
#define  MAX_RESTR_LEN                  7    /* max length of restr enzyme */
#define  DEFAULT_RESTR_FILE "restr_enzymes"
#define  MAX_PLIST                   10000    /* max # of sites for 1 enzyme*/

int    circular = 0;                          /* -c circular (loop) dna */
int    fragments = 0;                           /* write fragments (-f) */
int    sorted = 0;
int    nice_format = 0;
int    include_ends = 0;                      /* -e (for fragments, when linear */

char   filename[MAX_FILE_LEN];
char   restr_file[MAX_FILE_LEN];

int    pos_list[MAX_PLIST];                   /* list of re positions found */
int    pos_count;

static char restr_rcsid[] = 
      "$Id$";

/* integer comparison function for qsort() */

int int_compare( const void *a, const void *b )
   {
    if ( *(int *)a > *(int *)b)  return( 1 );
    if ( *(int *)a < *(int *)b)  return( -1 );
    return( 0 );
   }

void  print_list( int *list, int count )
   {
    int i;

    if ( nice_format )
       {
        printf( "\n" );
        for ( i = 0; i < count; i++ )
            printf( "%7d%c", list[i], ((i+1) % 6) ? ' ' : '\n' );
        printf( "\n" );
       }
    else
       {
        for ( i = 0; i < count; i++ )
            printf( "%7d", list[i] );
        printf( "\n" );
       }
   }

/* instead of printing positions, print_fragments prints the differences  */
/* between the positions. */

void  print_fragments( int len )
   {
    static int frags[MAX_PLIST+2];
    int        i, j;
    int        pos;

#ifdef DEBUG
    printf( "pf: pos_count %d\n", pos_count );
#endif
    if ( include_ends )
       {
        j = 0;
        pos = 0;
       }
    else
       {
        j = 1;
        pos = pos_list[0];
       }
        
    for ( i = 0; j < pos_count;  pos = pos_list[j++]  )
        frags[i++] =  pos_list[j] - pos;

    if ( include_ends )
        frags[i++] = len - pos;

    if ( sorted )
        qsort( (int *) frags, i, sizeof(int), int_compare );

    print_list( frags, i );
   }

/* add data from leaves to end of global array pos_list, updating pos_count  */
/* passed to and invoked from DFSprocess() for each subtrie corresponding to */
/* a restriction "hit"                                                       */

void  add_pos( TR_LEAF *l )
   {
    if ( l == NULL )
        return;
    pos_list[pos_count++] = (int) l->data;
    add_pos( l->next );
   }

void print_pos( int l )
   {
    printf( "%d\n", l );
   }

void  parse_args( int argc, char **argv )
   {
    int          c;
    int          i;
    extern char *optarg;
    extern int   optind;

    while ( (c = getopt(argc, argv, "cefns") ) != -1 )
        switch (c) 
           {
            case 'c':
                       circular = 1;
                       break;
            case 'e':
                       include_ends = 1;
                       break;
            case 'f':
                       fragments = 1;
                       break;
            case 'n':
                       nice_format = 1;
                       break;
            case 's':
                       sorted = 1;
                       break;
            case '?': 
                       fprintf( stderr, "Goodbye, Mr. Anderson\n" );
                       exit( 0 );
           }
    if ( optind >= argc )
       {
        fprintf( stderr, "need a sequence file name\n" );
        exit( 1 );
       }
    strncpy( filename, argv[optind], MAX_FILE_LEN );
    strncpy( restr_file, DEFAULT_RESTR_FILE, MAX_FILE_LEN );
   }

int  get_seq( int seqf, char *seq, int *len )
   {
    int            ret;
    static char    hdr[MAX_HDR_LENGTH];

    ret = read_seq( seqf, MAX_HDR_LENGTH, hdr, MAX_SEQ_LEN, seq );
    if ( ret > 0 ) 
       {
        *len = strlen( seq );
        uppercase( seq, seq );
       }
    return( ret );
   }


void  find_restr( TR_NODE *t, char *restr_file, int len )
    {
     FILE        *f;
     static char  line[MAX_LINE_LENGTH];
     static char  seq[MAX_RESTR_LEN];            /* enzyme seq */
     static char  seq2[MAX_RESTR_LEN];           /* enzyme seq */
     int          npats;
     static char  id[MAX_RESTR_LEN];             /* enzyme id */
     TR_NODE     *hit;

     if ( !( f = fopen( restr_file, "r" ) ) )
        {
         perror( restr_file );
         exit( errno );
        }
     while ( fgets( line, MAX_LINE_LENGTH, f ) )
        {
         npats = sscanf( line, "%s %s %s", id, seq, seq2 ) - 1;
         printf( "%8s ", id );
         pos_count = 0;                          /* reset positions array */ 
#ifdef DEBUG
         printf( "find: seq [%s]\n", seq );
#endif
         if ( hit = search_trie( t, seq, eqdna ) )
            {
#ifdef DEBUG
             printf( "got a hit:\n" );
             print_trie( hit, 0, print_pos );
#endif
             DFSprocess( hit->children, add_pos );
            }
         if ( npats == 2 )
             if ( hit = search_trie( t, seq2, eqdna ) )
                 DFSprocess( hit->children, add_pos );

         qsort( (int *) pos_list, pos_count, sizeof(int), int_compare );

         if ( fragments )
             print_fragments( len );
         else
             print_list( pos_list, pos_count );
        }
     fclose( f );
    }


main( int argc, char **argv )
   {
    int            seqf;
    static char    seq[MAX_SEQ_LEN];
    static char    tmp[MAX_SEQ_LEN];
    static char    tail[MAX_RESTR_LEN+1];
    TR_NODE       *t;
    char          *p;
    char          *seq_end;
    int            seq_len;

    parse_args( argc, argv );

    if ( ( seqf = open_seqf( filename, "acgtnxACGTNX", STD_IGNORED ) < 0 ) )
       {
        fprintf( stderr, "%s: fatal; ", progname );
        perror( filename );
        exit( errno );
       }
    if ( get_seq( seqf, seq, &seq_len ) <= 0 )
       {
        fprintf( stderr, "%s: no sequence in file %s\n", progname, filename );
        exit( 1 );
       }

    seq_end = seq + seq_len;             /* save end for circular case */
    if ( circular )
       {
        strncpy( tail, seq, MAX_RESTR_LEN );
        tail[MAX_RESTR_LEN] = '\0';
        strcat( seq, tail );             /* check for lengths!!!!!! */
       }

    t = create_trie();
    for ( p = seq; p < seq_end; p++ )
       {
        strncpy( tmp, p, MAX_RESTR_LEN );
#ifdef DEBUG
        printf( "add: [%s]\n", tmp );
#endif
        add_trie_string( t, tmp, (void *)( p - seq + 1 ) );
       }
#ifdef DEBUG
    print_trie( t, 0, print_pos );
#endif
    find_restr( t, restr_file, seq_len );
   }





