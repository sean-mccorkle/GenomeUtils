/* File:         neds.c                                                      */
/* Programmer:   Sean R. McCorkle, Brookhaven National Laboratory            */
/* Language:     C                                                           */
/*                                                                           */
/* Description:  NEDS - Neanderthal Embedded Documentation System            */
/*                                                                           */
/*               My attempt to embed man-page documentation IN the software  */
/*               in such a way that it can be easily converted into a        */
/*               reasonable text output, man page, or html page directly     */
/*               from the program.                                           */
/*                                                                           */
/* Usage:        include "neds.h"                                            */
/*                                                                           */
/*               neds_doc_out( mode, doc_string )                            */
/*                                                                           */
/*               where mode is one of                                        */
/*                                                                           */
/*                   NEDS_TEXT  for normal text output                       */
/*                   NEDS_MAN   for output into nroff -man for unix man page */
/*                   NEDS_HTML  for html formatted output                    */
/*                                                                           */
/*               and doc_string is a single string (including newlines) of   */
/*               the output document.                                        */
/*                                                                           */
/*               Lines begining with                                         */
/*                   __ are section headers (ie. OPTIONS, DESCRIPTION, etc)  */
/*                   _. are list items                                       */
/*                                                                           */
/*               Other lines are taken as text.  Text mode output pretty much*/
/*               mimics the actual docstring.  To preserve alignment in text */
/*               mode, if the first two chars are blanks, then they're       */
/*               stripped off.                                               */
/*                                                                           */
/*              In the text, constructions like <n> or <str> are converted   */
/*              to boldface n or str for man output, and blue font n and str */
/*                                                                           */
/*             for html output.                                              */
/*                                                                           */
/* Example:                                                                  */
/*                                                                           */
/*         neds_doc_out( "man", "\                                           */
/*__Name      pp - sequence pretty-printer                                \  */
/*                                                                        \  */
/*__Usage     pp [options] [<seq> [<seq> ... ] ]                          \  */
/*                                                                        \  */
/*            where <seq> is a file that contains one or more fasta format\  */
/*            DNA or amino acid sequences. A single hyphen (-) can be used\  */
/*            as a file name to specify stdin.  If no files are given,    \  */
/*            stdin is scanned.                                           \  */
/*                                                                        \  */
/*__Options                                                               \  */
/*_.           -B <n>   begin at position n                               \  */
/*_.           -L <n>   print <n> characters                              \  */
/*_.           -h       print help then exit                              \  */
/*_.           -l <n>   print <n> characters per line (default 50)        \  */
/*_.           -n       don't print position numbers on left              \  */
/*_.           -s <n>   insert a space every <n> characters (default 10), \  */
/*                      if 0 or <0, no subgrouping is done                \  */
/*_.           -t <str> use <str> as header title                         \  */
/*_.           -v       print version then exit                           \  */
/*__Note                                                                  \  */
/*           -B, -L and -t are applied to each sequence (they really      \  */
/*           only make sense for input which is a single sequence         \  */
/*              " )                                                          */
/*                                                                           */
/*                                                                           */
/* $Id$ */
/*****************************************************************************/

#include  <stdio.h>
#include  <string.h>
#include  "neds.h"

#define NEDS_STATE_CLEAR 0
#define NEDS_STATE_ITEM 1
int     neds_doc_html_state = CLEAR;

int   all_blank( char *start, char *stop );
void  neds_doc_header( int mode, char *title, char *date, char *etc );
void  neds_doc_line( int mode, char *start, char *stop );
void  neds_doc_sec_header( int mode, char *start, char *stop );
void  neds_doc_item( int mode, char *start, char *stop );
void  neds_doc_par( int mode );
void  clear_html_state( void );

/***************************************************/
/* This is the one routine intended for public use */
/***************************************************/

void  neds_doc_out( int mode, char *s )
   {
    char *nl;

    neds_doc_header( mode, "pp", "today", "etc" );
    while ( nl = index( s, '\n' ) )
       {
        neds_doc_line( mode, s, nl );
        s = nl + 1;
       }
    if ( *s )
       {
        nl = s + strlen( s );
        neds_doc_line( mode, s, nl );
       }
   }

/*************************************************************************/
/*  Here begin the private routines - not intended for use outside this  */
/*  module                                                               */
/*************************************************************************/

void  neds_doc_header( int mode, char *title, char *date, char *etc )
   {
    switch ( mode )
       {
        NEDS_TEXT:  printf( "%s %s %s\n", title, date, etc );
                    break;
        NEDS_MAN:   printf( ".TH \"%s\" \"%s\" \"%s\"\n", title, date, etc );
                    break;
        NEDS_HTML:  printf( "<title>%s</title>\n<blockquote>\n", title );
                    break;
       }

   }

/* this accepts an input line, decides what type it is, and then */
/* distributes it to the appropriate handling routine            */

void  neds_doc_line( int mode, char *start, char *stop )
   {
    if ( start[0] == "_" && start[1] == "_" )
        neds_doc_sec_header( mode, start+2, stop );
    else if ( start[0] == "_" && start[1] == "." )
        neds_doc_item( mode, start+2, stop );
    else if ( all_blank( start, stop ) )
        neds_doc_para( mode );
    else
        neds_doc_text( mod, start, stop );
   }


void  neds_doc_sec_header( int mode, char *start, char *stop )
   {
    local title, rest, sp1, sp2

                   title := tab( many( &letters ) )
                   sp1 := tab( many( ' ' ) )
                   rest := tab( 0 )
                   doc_sec_header( mode, title, rest, sp1 ) 
    switch ( mode )
       {qq
        NEDS_TEXT:  qwrite( title, sp, rest )
        NEDS_MAN:  {
                 write( ".SH ", map( title, &lcase, &ucase ) )
                 doc_text( mode, rest )
                }
        NEDS_HTML: {
                 clear_html_state()
                 write( "</blockquote>" )
                 write( "<h4>", title, "</h4>\n<blockquote>" )
                 doc_text( mode, rest )
                }
   }

void  neds_doc_item( int mode, char *start, char *stop )
    local title, rest, sp1, sp2


                   sp1 := tab( many( ' ' ) )
                   title := tab( upto( ' ' ) )
                   sp2 := tab( many( ' ' ) )
                   if ( match( "<" ) )
                      then {
                            title ||:= " " || tab( upto( ">" ) + 1 )
                            sp2 := tab( many( ' ' ) )
                           }
                   rest := tab( 0 )

    local t

    t := doc_process( mode, title )
    case mode of
       {
        NEDS_TEXT: write( sp1, title, sp2, rest )
        NEDS_MAN:  {
                 write( ".IP \"", t, "\"" )
                 doc_text( mode, rest )
	        }
        NEDS_HTML: {
                 set_html_state( "item" )
                 write( "<dt><tt>", t, "</tt></dt><dd>" )
                 doc_text( mode, rest )
	        }
       }
end

procedure doc_text( mode, s )
    case mode of
       {
        NEDS_MAN|NEDS_HTML:   s ? { tab( many( ' ' ) )
                              write( doc_process( mode, tab( 0 ) ) )
                            }
        NEDS_TEXT:         s ? { tab( match( "  " ) )
                              write( doc_process( mode, tab( 0 ) ) )
                            }
        default: write( s )
       }
end



/********************/
/* Paragraph marker */
/********************/

void  neds_doc_par( int mode )
   {
    switch ( mode )
       {
        NEDS_HTML:  clear_neds_html_state(); 
                    printf( "<P>\n" );
                    break;
        NEDS_MAN:   printf( ".PP\n" );
                    break;
        NEDS_TEXT:  printf( "\n" );
                    break;
       }
   }


#
# Convert <n> type constructions to fancy fonts if not text output
# this is probably better written with regexps

procedure doc_process( mode, s )

    local t, a

    if ( mode == "text" )
      then return( s )

    s ? { 
         t := ""
         while a := tab( upto( '<' ) ) do
	    {
             t ||:= a
             move( 1 )
             if any( &letters )
               then {
                     if  b := tab( upto( '<>' ) )
                       then {
                              if  match( ">" )
                                then
                                   t ||:= doc_trans( mode, b )
                                else
                                   t ||:= "<" || b || ">"
                              move( 1 )
		            }
                       else
                             t ||:= "<" || tab( 0 )
                    }
               else 
                     t ||:= "<"
	    }
         return( t || tab( 0 ))
        }
end

procedure  doc_trans( mode, s )
   case mode of
      {
       NEDS_HTML:  return( "<font color=\"blue\">" || s || "</font>" )
       NEDS_MAN:   return( "\\fI" || s || "\\fR" )
       default: return( s )
      }
end

void  set_neds_doc_html_state( int state )
   {
    if ( neds_doc_html_state == state )
        return;
    if ( neds_doc_html_state == NEDS_STATE_CLEAR && state == NEDS_STATE_ITEM )
        printf( "<DL>\n" );
    else
       {
        fprintf( stderr, "wierd html state!\n"  );
        exit( 1 );
       }
    neds_doc_html_state := state;
   }

void  clear_neds_doc_html_state( void )
   {
    if ( neds_doc_html_state == NEDS_STATE_ITEM )
        printf( "</DL>\n" );
    neds_doc_html_state = NEDS_STATE_CLEAR;
   }    


/*********************************************************************/
/* returns 1 (true) if the characters between start and stop are all */
/* blanks (or tabs), 0 (false) otherwise                             */
/*********************************************************************/

int  all_blank( char *start, char *stop )
   {
    while ( (start < stop ) && (*start == ' ' || *start == '\t' ) )
        start++;
    return ( start >= stop );
   }
