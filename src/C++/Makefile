OPTIONS= -O3
BINDIR=/usr/local/seq/bin

lossc:    lossc.o
	c++ $(OPTIONS) -o $@ $@.c++

install:
	chmod a+rx,g+w lossc post_lossc
	chgrp developer lossc post_lossc
	cp lossc post_lossc $(BINDIR)

.SUFFIXES: .c++ .o

.c++.o:
	c++ $(OPTIONS) -c $*.c++


