#
# makefile for Perple_X 6.9.2+
#
# To compile the Perple_X programs with this file first edit the compiler
# variables so that they are consistent with the fortran installation on 
# your system. Then type the command
#
#                             make -f makefile
#
# where makefile is the name of this file, e.g., LINUX_makefile.
#
# JADC, April 14, 2021    
# 
##################### COMPILER VARIABLES ####################################  
#                   
#                  FC = the name of the local fortran 77 compatible compiler
FC ?= gfortran
COMP77 ?= $(FC)
LINK77 ?= $(FC)


FFLAGS = -O3

LDFLAGS_RUNTIME :=

STATIC ?= 0
LDFLAGS_RUNTIME :=

ifeq ($(STATIC),1)
  # mostly static runtime everywhere (where supported)
  LDFLAGS_RUNTIME += -static-libgfortran -static-libgcc

  # Windows: try for a single self-contained .exe
  ifeq ($(EXT),.exe)
    LDFLAGS_RUNTIME += -static
    LDFLAGS_RUNTIME += -Wl,-Bstatic -lwinpthread -Wl,-Bdynamic
  endif

  # Linux: try to make quadmath static without using -static-libquadmath
  UNAME_S := $(shell uname -s 2>/dev/null)
  ifeq ($(UNAME_S),Linux)
    LDFLAGS_RUNTIME += -Wl,-Bstatic -lquadmath -Wl,-Bdynamic
  endif

endif

endif

FLINK += $(LDFLAGS_RUNTIME)

MYOBJ = actcor build fluids ctransf frendly meemum convex pstable pspts psvdraw pssect pt2curv vertex werami MC_fit
all: $(MYOBJ)

LIBRARY = pscom.o pslib.o rlib.o tlib.o flib.o olib.o resub.o minime_blas.o blas2lib.o

LIBRARY_PS = cont_lib.o pscom.o pslib.o flib.o tlib.o

clean: 
	rm -f *.o $(MYOBJ)

###################### TARGETS FOR 692+ PROGRAMS ######################### 
#
#
MC_fit: MC_fit.o $(LIBRARY)
	$(COMP77) $(FFLAGS) $(FLINK) $@.o $(LIBRARY) -o $@$(EXT)

convex: convex.o $(LIBRARY)
	$(COMP77) $(FFLAGS) $(FLINK) $@.o $(LIBRARY) -o $@$(EXT)

meemum: meemum.o $(LIBRARY)
	$(COMP77) $(FFLAGS) $(FLINK) $@.o $(LIBRARY) -o $@$(EXT)

vertex: vertex.o $(LIBRARY)
	$(COMP77) $(FFLAGS) $(FLINK) $@.o $(LIBRARY) -o $@$(EXT)

werami: werami.o $(LIBRARY)
	$(COMP77) $(FFLAGS) $(FLINK) $@.o $(LIBRARY) -o $@$(EXT)

pssect: psect.o cont_lib.o $(LIBRARY)
	$(COMP77) $(FFLAGS) $(FLINK) psect.o cont_lib.o $(LIBRARY) -o $@$(EXT)

pstable: pstable.o $(LIBRARY_PS)
	$(COMP77) $(FFLAGS) $(FLINK) $@.o $(LIBRARY_PS) -o $@$(EXT)

pspts: pspts.o $(LIBRARY_PS) 
	$(COMP77) $(FFLAGS) $(FLINK) $@.o $(LIBRARY_PS) -o $@$(EXT)

psvdraw: psvdraw.o $(LIBRARY_PS)
	$(COMP77) $(FFLAGS) $(FLINK) $@.o $(LIBRARY_PS) -o $@$(EXT)

pt2curv: pt2curv.o $(LIBRARY_PS)
	$(COMP77) $(FFLAGS) $(FLINK) $@.o $(LIBRARY_PS) -o $@$(EXT)

actcor: actcor.o $(LIBRARY)
	$(COMP77) $(FFLAGS) $(FLINK) $@.o $(LIBRARY) -o $@$(EXT)

build: build.o $(LIBRARY)
	$(COMP77) $(FFLAGS) $(FLINK) $@.o $(LIBRARY) -o $@$(EXT)

fluids: fluids.o $(LIBRARY)
	$(COMP77) $(FFLAGS) $(FLINK) $@.o $(LIBRARY) -o $@$(EXT)

ctransf: ctransf.o $(LIBRARY)
	$(COMP77) $(FFLAGS) $(FLINK) $@.o $(LIBRARY) -o $@$(EXT)

DEW_2_ver: DEW_2_ver.o $(LIBRARY)
	$(COMP77) $(FFLAGS) $(FLINK) $@.o $(LIBRARY) -o $@$(EXT)

frendly: frendly.o $(LIBRARY)
	$(COMP77) $(FFLAGS) $(FLINK) $@.o $(LIBRARY) -o $@$(EXT)


#################################################################################
# Explicitly make objects to override the default compiler choice implemented
# on some machines.

actcor.o: actcor.f perplex_parameters.h
	$(COMP77) $(FFLAGS) -c actcor.f
blas2lib.o: blas2lib.f perplex_parameters.h
	$(COMP77) $(FFLAGS) -c blas2lib.f
build.o: build.f perplex_parameters.h
	$(COMP77) $(FFLAGS) -c build.f
convex.o: convex.f perplex_parameters.h
	$(COMP77) $(FFLAGS) -c convex.f
cont_lib.o: cont_lib.f perplex_parameters.h
	$(COMP77) $(FFLAGS) -c cont_lib.f
ctransf.o: ctransf.f perplex_parameters.h
	$(COMP77) $(FFLAGS) -c ctransf.f
DEW_2_ver.o: DEW_2_ver.f perplex_parameters.h
	$(COMP77) $(FFLAGS) -c DEW_2_ver.f
flib.o: flib.f perplex_parameters.h
	$(COMP77) $(FFLAGS) -c flib.f
fluids.o: fluids.f perplex_parameters.h
	$(COMP77) $(FFLAGS) -c fluids.f
frendly.o: frendly.f perplex_parameters.h
	$(COMP77) $(FFLAGS) -c frendly.f
MC_fit.o: MC_fit.f perplex_parameters.h
	$(COMP77) $(FFLAGS) -c MC_fit.f
meemum.o: meemum.f perplex_parameters.h
	$(COMP77) $(FFLAGS) -c meemum.f
minime_blas.o: minime_blas.f perplex_parameters.h
	$(COMP77) $(FFLAGS) -c minime_blas.f
olib.o: olib.f perplex_parameters.h
	$(COMP77) $(FFLAGS) -c olib.f
pscom.o: pscom.f perplex_parameters.h
	$(COMP77) $(FFLAGS) -c pscom.f
psect.o: psect.f perplex_parameters.h
	$(COMP77) $(FFLAGS) -c psect.f
pslib.o: pslib.f perplex_parameters.h
	$(COMP77) $(FFLAGS) -c pslib.f
pspts.o: pspts.f perplex_parameters.h
	$(COMP77) $(FFLAGS) -c pspts.f
psvdraw.o: psvdraw.f perplex_parameters.h
	$(COMP77) $(FFLAGS) -c psvdraw.f
pstable.o: pstable.f perplex_parameters.h
	$(COMP77) $(FFLAGS) -c pstable.f
pt2curv.o: pt2curv.f perplex_parameters.h
	$(COMP77) $(FFLAGS) -c pt2curv.f
resub.o: resub.f perplex_parameters.h
	$(COMP77) $(FFLAGS) -c resub.f
rlib.o: rlib.f perplex_parameters.h
	$(COMP77) $(FFLAGS) -c rlib.f
tlib.o: tlib.f perplex_parameters.h
	$(COMP77) $(FFLAGS) -c tlib.f
vertex.o: vertex.f perplex_parameters.h
	$(COMP77) $(FFLAGS) -c vertex.f
werami.o: werami.f perplex_parameters.h
	$(COMP77) $(FFLAGS) -c werami.f

.c.o:
	$(CC) $(CFLAGS) $(INCL) -c $<

.f.o:
	$(FC) $(FFLAGS) $(INCL) -c $<
