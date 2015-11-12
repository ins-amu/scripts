include ../Makefile.inc

SOURCES = interface.cc \
        modelwriter.cc modelloader.cc trianglemesh.cc \
        helpers.cc averageplane.cc micropatch.cc microdelaunay.cc \
        relocation.cc polygon2.cc triangulator.cc \
        simplification.cc meshcleanup.cc meshoptimize.cc \
        oversampling.cc meshslicing.cc relaxation.cc relaxanglesmooth.cc \
        subdivbase.cc subdivlinear.cc subdivloop.cc meshdecimation.cc \
        edgeflips.cc delaunayflips.cc resampling.cc incdelaunay.cc \
        vertexinfo.cc vertexref.cc patch3d.cc patch2d_fb.cc patch2d.cc \
		densityfield.cc featureedges.cc meshskeleton.cc patchcache.cc \
		pntriangle.cc meshstats.cc cvtstats.cc gslpolysolve.cc

OBJECTS = $(foreach file,$(SOURCES),$(subst .cc,.o,$(file)))
LIBRARY = libremesh.a
GCC_INCL = -I..

all: ${OBJECTS}
	ar rcs ${LIBRARY} ${OBJECTS}
	chmod a+x ${LIBRARY}

test: ${OBJECTS} test.cc test.o
	${COMPILER} -o test test.o ${OBJECTS} -lpthread

multicore:
	$(MAKE) -j${CORES} all

depend:
	${COMPILER} ${GCC_INCL} -MM ${SOURCES} > Makefile.dep

%.o: %.cc
	${COMPILER} -c -o $@ $< ${GCC_FLAGS} ${GCC_INCL} -Wno-deprecated

clean: FORCE
	rm -f ${OBJECTS} ${LIBRARY} test.o test

FORCE:

include Makefile.dep
