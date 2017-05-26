#
#   1) make all
#   2) make install-local
#        or
#      make install
#   
#
# Change these according to your system or tastes
#
# for make install
PUBBINDIR   = /usr/local/seq/bin
# for make install-local   
LOCBINDIR   = ../bin

CC          = cc
MAXCOPT     = -O3

# for linux
# CLIBS       = 
# RANLIB      =
# for MacOSX
CLIBS       = -lm      
RANLIB      = ranlib   


CSRCS      = intervals.c kmers.c nt.c prosearch.c
#             io.c lpa_align.c nqcut.c repeats.c restr.c seqdiff.c sequtils.c \
#             trie.c sagetags.c atags.c gsts2.c lossc.c fcomp.c sageh.c \
#             overlap.c tagsearch.c sa_search.c intervals.c \
#             k-mer-directory.c
INCLUDES   =
#INCLUDES   = seqlib.h trie.h
SRCS       = $(CSRCS) $(INCLUDES)  Makefile
BINS       = intervals kmers nt prosearch
#BINS       = nt kmers repeats restr seqdiff nqcut sagetags atags gsts2 lossc fcomp \
#             sageh overlap prosearch tagsearch sa_search intervals k-mer-directory

COPT = $(MAXCOPT)


install: all
	- chgrp developer $(BINS) 
	- chmod a+rx,ug+w $(BINS)
	- mkdir -p $(PUBBINDIR)
	cp $(BINS) $(PUBBINDIR)

install-local: all
	- mkdir -p $(LOCBINDIR)
	cp $(BINS) $(LOCBINDIR)

all:  $(BINS)

kmers: kmers.o
	$(CC) $(COPTS) $(CCFLAGS) -o $@ $@.o

nt: nt.o
	$(CC) $(COPTS) $(CCFLAGS) -o $@ $@.o

#k-mer-directory:  k-mer-directory.o libseq.a
#	$(CC) $(COPT) $(CCFLAGS) -o $@ $@.o -L. -lseq -lm $(CLIBS)

intervals: intervals.o
	$(CC) $(COPTS) $(CCFLAGS) -o $@ $@.o

#tagsearch: tagsearch.o
#	$(CC) $(COPTS) $(CCFLAGS) -o $@ $@.o

prosearch: prosearch.o
	$(CC) $(COPTS) $(CCFLAGS) -o $@ $@.o

#overlap: overlap.o
#	$(CC) $(COPTS) $(CCFLAGS) -o $@ $@.o

#sagetags:  sagetags.o 
#	$(CC) $(COPT) $(CCFLAGS) -o $@ sagetags.o -lm $(CLIBS)

#atags:  atags.o   
#	$(CC) $(COPT) $(CCFLAGS) -o $@ atags.o
#
#n-mers:  n-mers.o libseq.a
#	$(CC) $(COPT) $(CCFLAGS) -o $@ n-mers.o -L. -lseq -lm $(CLIBS)
#
#sa_search:  sa_search.o libseq.a
#	$(CC) $(COPT) $(CCFLAGS) -o $@ sa_search.o -L. -lseq -lm $(CLIBS)
#
#fcomp:  fcomp.o libseq.a
#	$(CC) $(COPT) $(CCFLAGS) -o $@ fcomp.o -L. -lseq -lm $(CLIBS)
#
#gsts2:  gsts2.o   
#	$(CC) $(COPT) $(CCFLAGS) -o $@ gsts2.o
#
#lossc:  lossc.o
#	$(CC) $(COPT) $(CCFLAGS) -o $@ lossc.o $(CLIBS)
#
#cmap:   cmap.o  libseq.a
#	$(CC) $(COPT) $(CCFLAGS) -o $@ cmap.o -L. -lseq  $(CLIBS)
#
#nqcut:   nqcut.o  libseq.a
#	$(CC) $(COPT) $(CCFLAGS) -o $@ nqcut.o  -L. -lseq  $(CLIBS)
#
#repeats:   repeats.o  trie.o libseq.a
#	$(CC) $(COPT) $(CCFLAGS) -o $@ repeats.o trie.o -L. -lseq  $(CLIBS)
#
#restr:   restr.o  trie.o libseq.a
#	$(CC) $(COPT) $(CCFLAGS) -o $@ restr.o trie.o -L. -lseq  $(CLIBS)
#
#seqdiff:   seqdiff.o  libseq.a
#	$(CC) $(COPT) $(CCFLAGS) -o $@ seqdiff.o -L. -lseq  $(CLIBS)
#
#sageh:   sageh.o libseq.a
#	$(CC) $(COPT) $(CCFLAGS) -o $@ sageh.o -L. -lseq $(CLIBS)
#
#sagetr:   sagetr.o  trie.o libseq.a
#	$(CC) $(COPT) $(CCFLAGS) -o $@ sagetr.o trie.o -L. -lseq  $(CLIBS)
#
#libseq.a:  io.o sequtils.o lpa_align.o
#	ar r $@ io.o sequtils.o lpa_align.o
#	$(RANLIB) $@

.c.o:
	$(CC) -c $(COPT) $(CCFLAGS) $*.c

#seqdiff.o:    seqdiff.c   seqlib.h  
#io.o:         io.c        seqlib.h  
#sequtils.o:   sequtils.c  seqlib.h
#lpa_align.o:  lpa_align.c seqlib.h
#repeats.o:    repeats.c   seqlib.h  trie.h
#restr.o:      restr.c     seqlib.h  trie.h
#trie.o:       trie.c      seqlib.h  trie.h

clean:
	- rm $(BINS)
	- rm *.o


