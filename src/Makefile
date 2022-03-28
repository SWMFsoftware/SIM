# Copyright (C) 2002 Regents of the University of Michigan, portions used with
# permission. For more information, see http://csem.engin.umich.edu/tools/swmf

SHELL =/bin/sh

include ../Makefile.def

SEARCHDIR = -I${SHAREDIR}

include ../Makefile.conf
-include Makefile.DEPEND
-include Makefile.RULES

# Print environment variables for debugging.
print-% : ; @echo $* = $($*)

SEARCHDIR = -I${SHAREDIR}

# =============================================================================
# Define OBJECTS and MODULES

MODULES = \
	SIM_main.o \
	SIM_methods.o \
	SIM_solve.o \
	SIM_output.o


DEPEND:
	@perl ${SCRIPTDIR}/depend.pl ${SEARCHDIR} ${MODULES}

# =============================================================================
#
# Library for framework and also for stand alone code
#

MY_LIB  = libSIM.a

LIB: DEPEND
	make ${MY_LIB}
	@echo
	@echo ${MY_LIB} has been brought up to date.
	@echo

${MY_LIB}:	${MF} ${MODULES}
	rm -f ${MY_LIB}
	${AR} ${MY_LIB} ${MODULES}

# Object files not included into the library
OBJECTS_EXE = main.o

# Other required libraries
LIBSHARE  = ${LIBDIR}/libSHARE.a
LIBTIMING = ${LIBDIR}/libTIMING.a

# Libraries should be compiled first, because modules are used in main.
${OBJECTS_EXE}: ${LIBSHARE} ${LIBTIMING} ${MY_LIB}

EXE = ../SIM.exe

SIM:
	$(MAKE) ${EXE}
	@echo ' '
	@echo Program SIM has been brought up to date.
	@echo ' '

${EXE}: ${OBJECTS_EXE}
	rm -rf Tmp_; mkdir Tmp_
	cd Tmp_; \
		ar -x ../${MY_LIB}; \
		ar -x ${LIBTIMING}; \
		ar -x ${LIBSHARE}
	${LINK.f90} -o ${EXE} ${OBJECTS_EXE} Tmp_/*.o ${LBLAS} ${Lflag1} ${HYPRELIB}
	
# =============================================================================

clean: cleanfiles
	rm -rf Tmp_

distclean: clean

# keep this line
