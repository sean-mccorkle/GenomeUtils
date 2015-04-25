# C-level makefile for BNL Genome sequence utilities
#
# $Id: Makefile,v 0.1 2000/04/20 20:03:40 mccorkle Exp $
#

include ../Config/config.mak
include ../header.mak

CSRCS      = seqdiff.c io.c lpa_align.c sequtils.c
INCLUDES   = seqlib.h
SRCS       = $(CSRCS) $(INCLUDES)  Makefile
BINS       = seqdiff

all:  $(BINS)


seqdiff:   seqdiff.o  libseq.a
	cc $(CCFLAGS) -o $@ seqdiff.o -L. -lseq  $(CLIBS)

sagetr:   sagetr.o  trie.o libseq.a
	cc $(CCFLAGS) -o $@ sagetr.o trie.o -L. -lseq  $(CLIBS)

hairpin:   hairpin.o  trie.o libseq.a
	cc $(CCFLAGS) -o $@ hairpin.o trie.o -L. -lseq  $(CLIBS)

libseq.a:  io.o sequtils.o lpa_align.o
	ar r $@ io.o sequtils.o lpa_align.o

.c.o:
	cc -c $(CCFLAGS) $*.c

seqdiff.o:    seqdiff.c   seqlib.h  
io.o:         io.c        seqlib.h  
sequtils.o:   sequtils.c  seqlib.h
lpa_align.o:  lpa_align.c seqlib.h




