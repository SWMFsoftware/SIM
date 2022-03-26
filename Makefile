# Copyright (C) 2002 Regents of the University of Michigan, portions used with
# permission. For more information, see http://csem.engin.umich.edu/tools/swmf

DEFAULT_TARGET = SIM
DEFAULT_EXE    = ${DEFAULT_TARGET}.exe

default : ${DEFAULT_TARGET}

include Makefile.def
include Makefile.conf

help:
	@echo 'Please write something here.'

info:
	@echo "Total lines of Fortran: `wc -l src*/*.f* share/Library/src/*.f* | tail -1`"

install:
	touch src/Makefile.DEPEND

LIB:
	cd src; $(MAKE) LIB
	cd srcInterface; $(MAKE) LIB

SIM:
	cd ${SHAREDIR}; $(MAKE) LIB
	cd ${TIMINGDIR}; $(MAKE) LIB
	cd src; $(MAKE) LIB
	cd src; make SIM


# Default component
COMPONENT = IE

rundir:
	mkdir -p ${RUNDIR}/IE/output

clean:
	cd src; make clean

distclean:
	./Config.pl -uninstall

allclean:
	cd src; make distclean
	rm -f *~

test:
	echo "There is no test for IE/SIMPLE yet." > notest.diff
