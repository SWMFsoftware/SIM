# Copyright (C) 2002 Regents of the University of Michigan, portions used with
# permission. For more information, see http://csem.engin.umich.edu/tools/swmf

SHELL =/bin/sh

include ../Makefile.def
include ../Makefile.conf
-include Makefile.DEPEND

# Configurable objects

OBJECTS = \
	IE_wrapper.o

MODULES = \
	IE_ModUnused.o

SEARCHDIR = -I${SHAREDIR} -I${COUPLERDIR} -I${LIBRARYDIR} -I../src

DEPEND:
	@perl ${SCRIPTDIR}/depend.pl ${SEARCHDIR} ${MODULES} ${OBJECTS}

#
#	Making executables and libraries
#
MY_LIB = ${LIBDIR}/libIE.a

LIB:	DEPEND
	make ${MY_LIB}
	@echo
	@echo ${MY_LIB} has been brought up to date.
	@echo

SIM_LIB = ../src/libSIM.a 

${MY_LIB}: ${SIM_LIB} ${MODULES} ${OBJECTS}
	cp -f ${SIM_LIB} ${MY_LIB}
	${AR} ${MY_LIB} ${MODULES} ${OBJECTS}

clean: cleanfiles

distclean: clean
	rm -f Makefile.DEPEND

