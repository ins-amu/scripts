include ../Makefile.inc

SOURCES = main.cc
OBJECTS = $(foreach file,$(SOURCES),$(subst .cc,.o,$(file)))
BINARY = cmdremesher

GCC_INCL = -I..
GCC_LINK = -L../libremesh -lremesh -lpthread
GSL_LIBS = $(shell pkg-config --libs gsl)

all: libraries ${OBJECTS}
	${COMPILER} -o ${BINARY} ${OBJECTS} ${GCC_LINK} ${GSL_LIBS}

multicore:
	$(MAKE) -j${CORES} all

depend:
	${COMPILER} ${GCC_INCL} -MM ${SOURCES} > Makefile.dep

%.o: %.cc
	${COMPILER} -c -o $@ $< ${GCC_FLAGS} ${GCC_INCL}

clean: FORCE
	rm -f ${BINARY} ${OBJECTS}

libraries: FORCE
	$(MAKE) -C ../libremesh

FORCE:

include Makefile.dep
