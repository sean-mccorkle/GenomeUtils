#
# Common code for makefiles for BNL Genome sequence utilties package
#

install: all
	- chgrp developer $(BINS) 
	- chmod a+rx,ug+w $(BINS) 
	cp $(BINS) $(PUBBINDIR)

install-local: all
	cp $(BINS) $(LOCBINDIR)

dist:
	echo making  dist, dir is $(DIR)
	cd ..; tar rvf ../$(DISTFILE) \
	    `echo $(SRCS) | tr ' ' '\n' | sed -e 's/^/$(DIR)\//'`

clean:
	rm -f *~ #*#

showrcs:
	for p in $(SRCS) ; \
	    do rcsdiff $$p ; \
	    done


