/* File:         neds.c                                                      */
/* Programmer:   Sean R. McCorkle, Brookhaven National Laboratory            */
/* Language:     C                                                           */
/*                                                                           */
/* Description:  Header file for NEDS:  Neanderthal Embedded Documentation   */
/*               System (see neds.c for more info)                           */
/*                                                                           */
/* $Id$ */
/*****************************************************************************/


void neds_doc_out ( int mode, char *doc_string );

/* Allowed values for output mode */

#define  NEDS_TEXT  1   /* for normal text output                        */
#define  NEDS_MAN   2   /*  for output into nroff -man for unix man page */
#define  NEDS_HTML  3   /* for html formatted output                     */


