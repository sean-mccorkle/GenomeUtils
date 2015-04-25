# C-level makefile for BNL Genome sequence utilities
#
# $Id: Makefile,v 0.4 2003/11/10 21:01:30 mccorkle Exp mccorkle $
#

include ../Config/config.mak
include ../header.mak

CSRCS      = io.c lpa_align.c nqcut.c repeats.c restr.c seqdiff.c sequtils.c \
             trie.c sagetags.c atags.c gsts2.c lossc.c fcomp.c sageh.c \
             overlap.c prosearch.c tagsearch.c sa_search.c intervals.c
INCLUDES   = seqlib.h trie.h
SRCS       = $(CSRCS) $(INCLUDES)  Makefile
BINS       = repeats restr seqdiff nqcut sagetags atags gsts2 lossc fcomp \
             sageh overlap prosearch tagsearch sa_search intervals

COPT = $(MAXCOPT)

all:  $(BINS)

intervals: intervals.o
	$(CC) $(COPTS) $(CCFLAGS) -o $@ $@.o

tagsearch: tagsearch.o
	$(CC) $(COPTS) $(CCFLAGS) -o $@ $@.o

prosearch: prosearch.o
	$(CC) $(COPTS) $(CCFLAGS) -o $@ $@.o

overlap: overlap.o
	$(CC) $(COPTS) $(CCFLAGS) -o $@ $@.o

sagetags:  sagetags.o 
	$(CC) $(COPT) $(CCFLAGS) -o $@ sagetags.o -lm $(CLIBS)

atags:  atags.o   
	$(CC) $(COPT) $(CCFLAGS) -o $@ atags.o

n-mers:  n-mers.o libseq.a
	$(CC) $(COPT) $(CCFLAGS) -o $@ n-mers.o -L. -lseq -lm $(CLIBS)

sa_search:  sa_search.o libseq.a
	$(CC) $(COPT) $(CCFLAGS) -o $@ sa_search.o -L. -lseq -lm $(CLIBS)

fcomp:  fcomp.o libseq.a
	$(CC) $(COPT) $(CCFLAGS) -o $@ fcomp.o -L. -lseq -lm $(CLIBS)

gsts2:  gsts2.o   
	$(CC) $(COPT) $(CCFLAGS) -o $@ gsts2.o

lossc:  lossc.o
	$(CC) $(COPT) $(CCFLAGS) -o $@ lossc.o $(CLIBS)

cmap:   cmap.o  libseq.a
	$(CC) $(COPT) $(CCFLAGS) -o $@ cmap.o -L. -lseq  $(CLIBS)

nqcut:   nqcut.o  libseq.a
	$(CC) $(COPT) $(CCFLAGS) -o $@ nqcut.o  -L. -lseq  $(CLIBS)

repeats:   repeats.o  trie.o libseq.a
	$(CC) $(COPT) $(CCFLAGS) -o $@ repeats.o trie.o -L. -lseq  $(CLIBS)

restr:   restr.o  trie.o libseq.a
	$(CC) $(COPT) $(CCFLAGS) -o $@ restr.o trie.o -L. -lseq  $(CLIBS)

seqdiff:   seqdiff.o  libseq.a
	$(CC) $(COPT) $(CCFLAGS) -o $@ seqdiff.o -L. -lseq  $(CLIBS)

sageh:   sageh.o libseq.a
	$(CC) $(COPT) $(CCFLAGS) -o $@ sageh.o -L. -lseq $(CLIBS)

sagetr:   sagetr.o  trie.o libseq.a
	$(CC) $(COPT) $(CCFLAGS) -o $@ sagetr.o trie.o -L. -lseq  $(CLIBS)

libseq.a:  io.o sequtils.o lpa_align.o
	ar r $@ io.o sequtils.o lpa_align.o
	$(RANLIB) $@

.c.o:
	$(CC) -c $(COPT) $(CCFLAGS) $*.c

seqdiff.o:    seqdiff.c   seqlib.h  
io.o:         io.c        seqlib.h  
sequtils.o:   sequtils.c  seqlib.h
lpa_align.o:  lpa_align.c seqlib.h
repeats.o:    repeats.c   seqlib.h  trie.h
restr.o:      restr.c     seqlib.h  trie.h
trie.o:       trie.c      seqlib.h  trie.h


