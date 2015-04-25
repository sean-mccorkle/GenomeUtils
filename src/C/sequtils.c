/* Module:       sequtils.c                                               */
/* Programmer:   Sean R. McCorkle                                         */
/* Lanuage:      C                                                        */
/* Description:  Small, self-contained routines for manipulating sequence */
/*               strings                                                  */
/*                                                                        */
/* $Id: sequtils.c,v 3.5 2007/11/13 19:30:28 mccorkle Exp mccorkle $      */
/**************************************************************************/


static char sequtils_rcs_id[] =
      "$Id: sequtils.c,v 3.5 2007/11/13 19:30:28 mccorkle Exp mccorkle $";

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "seqlib.h"

                              /***********/
                              /* Globals */
                              /***********/


char  *progname = "";                        /* for error messages        */

/* Used by char_to_sequential to quickly map ascii characters to a 1-15 */
/* sequential code */


char seq_map[] = {
                  /*          0     1     2     3     4     5     6     7   */
                  /* 000 */   0,    X,    X,    X,    X,    X,    X,    X,
                  /* 010 */   X,    X,    X,    X,    X,    X,    X,    X,
                  /* 020 */   X,    X,    X,    X,    X,    X,    X,    X,
                  /* 030 */   X,    X,    X,    X,    X,    X,    X,    X,
                  /* 040 */   X,    X,    X,    X,    X,    X,    X,    X,
                  /* 050 */   X,    X,    X,    X,    X,    X,    X,    X,
                  /* 060 */   X,    X,    X,    X,    X,    X,    X,    X,
                  /* 070 */   X,    X,    X,    X,    X,    X,    X,    X,
                  /* 100 */   X,    1,   14,    2,   13,    X,    X,    3,
                  /* 110 */  12,    X,    X,   10,    X,    5,   15,    X,
                  /* 120 */   X,    X,    6,    8,    4,    X,   11,    7,
                  /* 130 */   0,    9,    X,    X,    X,    X,    X,    X,
                  /* 140 */   X,    1,   14,    2,   13,    X,    X,    3,
                  /* 150 */  12,    X,    X,   10,    X,    5,   15,    X,
                  /* 160 */   X,    X,    6,    8,    4,    X,   11,    7,
                  /* 170 */   0,    9,    X,    X,    X,    X,    X,    X,
                 };

char agree[16][16] = {
/*         a   c   g   t   m   r   w   s   y   k   v   h   d   b   n */
       0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
/*a*/  0,  1,  0,  0,  0,  1,  1,  1,  0,  0,  0,  1,  1,  1,  0,  1,
/*c*/  0,  0,  1,  0,  0,  1,  0,  0,  1,  1,  0,  1,  1,  0,  1,  1,
/*g*/  0,  0,  0,  1,  0,  0,  1,  0,  1,  0,  1,  1,  0,  1,  1,  1,
/*t*/  0,  0,  0,  0,  1,  0,  0,  1,  0,  1,  1,  0,  1,  1,  1,  1,
/*m*/  0,  1,  1,  0,  0,  1,  1,  1,  1,  1,  0,  1,  1,  1,  1,  1,
/*r*/  0,  1,  0,  1,  0,  1,  1,  1,  1,  0,  1,  1,  1,  1,  1,  1,
/*w*/  0,  1,  0,  0,  1,  1,  1,  1,  0,  1,  1,  1,  1,  1,  1,  1,
/*s*/  0,  0,  1,  1,  0,  1,  1,  0,  1,  1,  1,  1,  1,  1,  1,  1,
/*y*/  0,  0,  1,  0,  1,  1,  0,  1,  1,  1,  1,  1,  1,  1,  1,  1,
/*k*/  0,  0,  0,  1,  1,  0,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,
/*v*/  0,  1,  1,  1,  0,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,
/*h*/  0,  1,  1,  0,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,
/*d*/  0,  1,  0,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,
/*b*/  0,  0,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,
/*n*/  0,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,
};


/* char_to_sequential() translates a DNA character sequence, character-by-  */
/* character, to a sequential mapping, according to the scheme              */
/* a -> 1, c-> 2, g -> 3, t -> 4, m -> 5, r -> 6, w -> 7, s -> 8, y -> 9,   */
/* k -> 10, v -> 11, h -> 12, d -> 13, b -> 14, n -> 15                     */
/* Note: the calling routine is responsible for string allocation           */

void  char_to_sequential( char *characters, char *sequential )
   {
    while ( *characters )
       {
        *sequential++ = SEQUENT( *characters ); 
        characters++;
       }
    *sequential = '\0';
   }

/* sequential_to_char() does the reverse mapping. Note: only returns  */
/* lowercase                                                          */

void  sequential_to_char( char *sequential, char *characters )
   {
    while ( *sequential )
       {
        *characters++ = CHARAC( *sequential );
        sequential++;
       }
    *characters = '\0';
   }

/* comparision functions (useful for trie) */

int  eqchar( char a, char b )   
   {
    return ( a == b );
   }

int  eqdna( char a, char b )
   {
    return( AGREE( SEQUENT(a), SEQUENT(b) )  );
   }


/* Used by dna_reverse_comp to quickly map nucleotide characters and     */
/* ambiguity codes to the appropriate complement                         */

char comp_char[] = {
                  /*          0     1     2     3     4     5     6     7   */
                  /* 000 */   X,    X,    X,    X,    X,    X,    X,    X,
                  /* 010 */   X,    X,    X,    X,    X,    X,    X,    X,
                  /* 020 */   X,    X,    X,    X,    X,    X,    X,    X,
                  /* 030 */   X,    X,    X,    X,    X,    X,    X,    X,
                  /* 040 */   X,    X,    X,    X,    X,    X,    X,    X,
                  /* 050 */   X,    X,    X,    X,    X,    X,    X,    X,
                  /* 060 */   X,    X,    X,    X,    X,    X,    X,    X,
                  /* 070 */   X,    X,    X,    X,    X,    X,    X,    X,
                  /* 100 */   X,   'T',  'V',  'G',  'H',   X,    X,   'C',
                  /* 110 */  'D',   X,    X,   'M',   X,   'K',  'N',   X,
                  /* 120 */   X,    X,   'Y',  'S',  'A',   X,   'B',  'W',
                  /* 130 */   X,   'R',   X,    X,    X,    X,    X,    X,
                  /* 140 */   X,   't',  'v',  'g',  'h',   X,    X,   'c',
                  /* 150 */  'd',   X,    X,   'm',   X,   'k',  'n',   X,
                  /* 160 */   X,    X,   'y',  's',  'a',   X,   'b',  'w',
                  /* 170 */   X,   'r',   X,    X,    X,    X,    X,    X,
                 };


void  dna_reverse_comp( char *s, char *r )
   {
    int  l;
    char *t;

    l = strlen( s );
    t = r + l;
    *t-- = '\0';
    while ( *s )
        *t-- = comp_char[ *s++ ];
   }



/* seq_reverse_comp() - like dna_reverse_comp() except it works on sequential*/
/* mapped strings.                                                           */

char seq_comp[] = { X, 4, 3, 2, 1, 10, 9, 7, 8, 6, 5, 14, 13, 12, 11, 15 };

void  seq_reverse_comp( char *s, char *r )
   {
    int  l;
    char *t;

    l = strlen( s );
    t = r + l;
    *t-- = '\0';
    while ( *s )
        *t-- = seq_comp[ *s++ ];
   }

/* blank() returns 1 if the input string is all whitespace, 0 otherwise */

int blank( char *s )

   {
    while ( *s )
       {
        if ( *s != ' ' && *s != '\t' && *s != '\n' )
            return( 0 );
        s++;
       }
    return( 1 );
   }


static struct trans_tab_t { char *codon; char aa; }
    trans_tab[] = { 
       { "AAA", 'K' }, { "AAC", 'N' }, { "AAG", 'K' }, { "AAT", 'N' },
       { "ACA", 'T' }, { "ACC", 'T' }, { "ACG", 'T' }, { "ACT", 'T' },
       { "AGA", 'R' }, { "AGC", 'S' }, { "AGG", 'R' }, { "AGT", 'S' },
       { "ATA", 'I' }, { "ATC", 'I' }, { "ATG", 'M' }, { "ATT", 'I' },

       { "CAA", 'Q' }, { "CAC", 'H' }, { "CAG", 'Q' }, { "CAT", 'H' },
       { "CCA", 'P' }, { "CCC", 'P' }, { "CCG", 'P' }, { "CCT", 'P' },
       { "CGA", 'R' }, { "CGC", 'R' }, { "CGG", 'R' }, { "CGT", 'R' },
       { "CTA", 'L' }, { "CTC", 'L' }, { "CTG", 'L' }, { "CTT", 'L' },
                    
       { "GAA", 'E' }, { "GAC", 'D' }, { "GAG", 'E' }, { "GAT", 'D' },
       { "GCA", 'A' }, { "GCC", 'A' }, { "GCG", 'A' }, { "GCT", 'A' },
       { "GGA", 'G' }, { "GGC", 'G' }, { "GGG", 'G' }, { "GGT", 'G' },
       { "GTA", 'V' }, { "GTC", 'V' }, { "GTG", 'V' }, { "GTT", 'V' },

       { "TAA", '*' }, { "TAC", 'Y' }, { "TAG", '*' }, { "TAT", 'Y' },
       { "TCA", 'S' }, { "TCC", 'S' }, { "TCG", 'S' }, { "TCT", 'S' },
       { "TGA", '*' }, { "TGC", 'C' }, { "TGG", 'W' }, { "TGT", 'C' },
       { "TTA", 'L' }, { "TTC", 'F' }, { "TTG", 'L' }, { "TTT", 'F' }
     };


char trans_codon( char n1, char n2, char n3 )
   {
    int   i;
    static char cod[4];
    cod[0] = n1;
    cod[1] = n2;
    cod[2] = n3;
    cod[3] = '\0';

    i = 0;
    while ( i < 64 && strcmp( cod, trans_tab[i].codon ) != 0 )
       i++;
    return(  i < 64 ? trans_tab[i].aa : '?' );
   }

/**************************************************************************/
/* puts protein translation of dna[] into aa[] in triplicate              */
/* form, so that the translated string has the same length as the DNA     */
/* string, and each base maps to its amino-acid at the same position.     */
/*  i.e.    "GATCCAGCG"  translates to "DDDPPPAAA", not "DPA".            */
/* translation is assumed to start in reading frame 1.  Any leftover      */
/* bases at the end translate to "?"                                      */
/* aa[] must already be allocated and must be the same size as dna[]      */
/**************************************************************************/

void  translate( char *dna, char *aa )
   {
    int  len;
    int  l;
    int  i;
    char a;
    len = strlen( dna );
    l = (len / 3) * 3;   
    for ( i = 0; i < l; i += 3 )
       {
        a = trans_codon( dna[i], dna[i+1], dna[i+2] );
        aa[i] = a;
        aa[i+1] = a;
        aa[i+2] = a;
       }
    while ( i < len )
       aa[i++] = '?';
    aa[i] = '\0';
   }



/*****************************************************************/
/* lowercase() copies characters from the array u to the array   */
/* l, converting them to lower case as it does it.  The calling  */
/* routine is responsible for array allocation.                  */
/*****************************************************************/

void  lowercase( char *u, char *l )

   {
    while( *u )
       {
        *l++ = tolower( *u );
        u++;
       }
    *l = '\0';
   }

/*****************************************************************/
/* uppercase() copies characters from the array l to the array   */
/* u, converting them to upper case as it does it.  The calling  */
/* routine is responsible for array allocation.                  */
/*****************************************************************/

void  uppercase( char *l, char *u )

   {
    while( *l ) 
       {
        *u++ = toupper( *l );
        l++;
       }
    *u = '\0';
   }


/***********************************************************************/
/* This function effectively removes any directory paths and returns a */
/* pointer to the last file part of a full path name.                  */
/***********************************************************************/

char *remove_path( char *path )

   {
    char *s;

    for ( s = path; *s; s++ )               /* move s to the end '\0'        */
       ;
    while ( s < path && *(s-1) != '/' )     /* backup until '/' or beginning */
        s--;
    return( s );                            /* voila                         */
   }





/* Hmm.  Seems like MetroWerks C ansi library is missing strdup().  */

#ifdef METROWERKS
char *strdup( char *s )
   {
    char *t;

    t = (char *) malloc_safely( strlen( s ) + 1, "occurred in strdup" );
    strcpy( t, s );
    return( t );
   }
#endif


