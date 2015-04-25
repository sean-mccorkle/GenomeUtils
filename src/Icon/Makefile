# Makefile for BNL Genome Icon package
#
# $Id: Makefile,v 0.6 2000/04/26 15:52:29 mccorkle Exp mccorkle $
#

include ../Config/config.mak
include ../header.mak

BINS      = orfs orfs2tv pp ranseq res sage sre trans tview
SRCS      = orfs.icn orfs2tv.icn pp.icn ranseq.icn res.icn \
            sage.icn sre.icn trans.icn tview.icn \
            graf.icn io.icn doc.icn \
            Makefile


all:  $(BINS)


orfs:  orfs.icn io.u1 io.u2
	icont $@

orfs2tv:  orfs2tv.icn io.u1 io.u2
	icont $@

pp:      pp.icn io.u1 io.u2 doc.u1 doc.u2
	icont $@

ranseq:  ranseq.icn io.u1  io.u2
	icont $@

rc:      rc.icn io.u1  io.u2
	icont $@

res:     res.icn io.u1  io.u2
	icont $@

sage:   sage.icn io.u1  io.u2
	icont $@

s2:   s2.icn io.u1  io.u2
	icont $@

sre:  sre.icn io.u1 io.u2
	icont $@

splt:   splt.icn io.u1  io.u2
	icont $@

trans:  trans.icn io.u1 io.u2
	icont $@

tview:  tview.icn  graf.u1  graf.u2
	icont $@

#iceu:   iceu.icn
#	icont $@

.SUFFIXES :  .icn .u1 .u2

.icn.u1:
	icont -c $<
.icn.u2:
	icont -c $<


