# Top level makefile for BNL Genome sequence utilities packages
#
# $Id: Makefile,v 0.6 2000/03/21 20:38:57 mccorkle Exp mccorkle $
#

DISTFILE    = gensrcdist.tar
SUBDIRS     = Scripts Icon C
OSLIST      = solaris|irix6|linux|openbsd|macosx|pando


all install install-local showrcs:
	for d in $(SUBDIRS)                                             ; \
	    do ( cd $$d; make $@ DISTFILE=$(DISTFILE) )                 ; \
            done

conf:
	( cd Config; make conf OS=$(OS) )

clean:
	rm -f *~ #*#
	for d in $(SUBDIRS)                                             ; \
	    do ( cd $$d; make $@ )                                      ; \
            done

dist:
	tar cvf ../$(DISTFILE) Config Makefile *.mak   # documents too?
	for d in $(SUBDIRS)                                             ; \
	    do ( cd $$d; make dist DIR=$$d DISTFILE=$(DISTFILE) )       ; \
            done

install-remote: dist
	cat Config/machines
	for m in `cat Config/machines`                                  ; \
	    do                                                            \
            case $$m in                                                   \
                $(OSLIST) ) os=$$m ;;                                     \
	        *)  echo distributing to $$m, an $$os machine           ; \
                    ssh $$m mkdir /tmp/gensrc                           ; \
	            scp ../$(DISTFILE) $$m:/tmp/gensrc                  ; \
	            ssh $$m "cd /tmp/gensrc;                              \
                             tar xvf $(DISTFILE);                         \
                             make conf OS=$$os;                         \
                             make install"                              ; \
	            ssh $$m rm -r /tmp/gensrc  ;;                         \
              esac                                                      ; \
            done
	echo remote distribution completed






