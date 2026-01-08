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
STATIC ?= 0

# Default executable extension
EXT ?=
ifeq ($(OS),Windows_NT)
  EXT := .exe
endif

FFLAGS = -O3

LDFLAGS_RUNTIME_PRE  :=
LDFLAGS_RUNTIME_POST :=

ifeq ($(STATIC),1)
  # remove libgfortran + libgcc* DLL/SO deps
  LDFLAGS_RUNTIME_PRE += -static-libgfortran -static-libgcc

  # key for Linux: prevents gfortran-injected shared libs (like -lquadmath)
  # from being recorded as NEEDED when they are not actually required
  LDFLAGS_RUNTIME_PRE += -Wl,--as-needed

  UNAME_S := $(shell uname -s 2>/dev/null)

  # Windows (MSYS2/MinGW/Cygwin): kill libquadmath-0.dll and libwinpthread-1.dll
  ifneq (,$(filter MINGW% MSYS% CYGWIN%,$(UNAME_S)))
    QUADMATH_A := $(shell $(FC) -print-file-name=libquadmath.a)
    WINPTHREAD_A := $(shell $(FC) -print-file-name=libwinpthread.a)

    ifneq ($(QUADMATH_A),libquadmath.a)
      LDFLAGS_RUNTIME_POST += -Wl,-Bstatic $(QUADMATH_A) -Wl,-Bdynamic
    else
      LDFLAGS_RUNTIME_POST += -Wl,-Bstatic -lquadmath -Wl,-Bdynamic
    endif

    ifneq ($(WINPTHREAD_A),libwinpthread.a)
      LDFLAGS_RUNTIME_POST += -Wl,-Bstatic $(WINPTHREAD_A) -Wl,-Bdynamic
    else
      LDFLAGS_RUNTIME_POST += -Wl,-Bstatic -lwinpthread -Wl,-Bdynamic
    endif
  endif

  # Linux: kill libquadmath.so.0
  ifeq ($(UNAME_S),Linux)
    QUADMATH_A := $(shell $(FC) -print-file-name=libquadmath.a)
    ifneq ($(QUADMATH_A),libquadmath.a)
      LDFLAGS_RUNTIME_POST += -Wl,-Bstatic $(QUADMATH_A) -Wl,-Bdynamic
    else
      LDFLAGS_RUNTIME_POST += -Wl,-Bstatic -lquadmath -Wl,-Bdynamic
    endif
  endif
endif

FLINK_PRE  += $(LDFLAGS_RUNTIME_PRE)
FLINK_POST += $(LDFLAGS_RUNTIME_POST)

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
	$(COMP77) $(FFLAGS) $(FLINK_PRE) $@.o $(LIBRARY) -o $@$(EXT) $(FLINK_POST)

convex: convex.o $(LIBRARY)
	$(COMP77) $(FFLAGS) $(FLINK_PRE) $@.o $(LIBRARY) -o $@$(EXT) $(FLINK_POST)

meemum: meemum.o $(LIBRARY)
	$(COMP77) $(FFLAGS) $(FLINK_PRE) $@.o $(LIBRARY) -o $@$(EXT) $(FLINK_POST)

vertex: vertex.o $(LIBRARY)
	$(COMP77) $(FFLAGS) $(FLINK_PRE) $@.o $(LIBRARY) -o $@$(EXT) $(FLINK_POST)

werami: werami.o $(LIBRARY)
	$(COMP77) $(FFLAGS) $(FLINK_PRE) $@.o $(LIBRARY) -o $@$(EXT) $(FLINK_POST)

pssect: psect.o cont_lib.o $(LIBRARY)
	$(COMP77) $(FFLAGS) $(FLINK) psect.o cont_lib.o $(LIBRARY) -o $@$(EXT)

pstable: pstable.o $(LIBRARY_PS)
	$(COMP77) $(FFLAGS) $(FLINK_PRE) $@.o $(LIBRARY_PS) -o $@$(EXT) $(FLINK_POST)

pspts: pspts.o $(LIBRARY_PS) 
	$(COMP77) $(FFLAGS) $(FLINK_PRE) $@.o $(LIBRARY_PS) -o $@$(EXT) $(FLINK_POST)

psvdraw: psvdraw.o $(LIBRARY_PS)
	$(COMP77) $(FFLAGS) $(FLINK_PRE) $@.o $(LIBRARY_PS) -o $@$(EXT) $(FLINK_POST)

pt2curv: pt2curv.o $(LIBRARY_PS)
	$(COMP77) $(FFLAGS) $(FLINK_PRE) $@.o $(LIBRARY_PS) -o $@$(EXT) $(FLINK_POST)

actcor: actcor.o $(LIBRARY)
	$(COMP77) $(FFLAGS) $(FLINK_PRE) $@.o $(LIBRARY) -o $@$(EXT) $(FLINK_POST)

build: build.o $(LIBRARY)
	$(COMP77) $(FFLAGS) $(FLINK_PRE) $@.o $(LIBRARY) -o $@$(EXT) $(FLINK_POST)

fluids: fluids.o $(LIBRARY)
	$(COMP77) $(FFLAGS) $(FLINK_PRE) $@.o $(LIBRARY) -o $@$(EXT) $(FLINK_POST)

ctransf: ctransf.o $(LIBRARY)
	$(COMP77) $(FFLAGS) $(FLINK_PRE) $@.o $(LIBRARY) -o $@$(EXT) $(FLINK_POST)

DEW_2_ver: DEW_2_ver.o $(LIBRARY)
	$(COMP77) $(FFLAGS) $(FLINK_PRE) $@.o $(LIBRARY) -o $@$(EXT) $(FLINK_POST)

frendly: frendly.o $(LIBRARY)
	$(COMP77) $(FFLAGS) $(FLINK_PRE) $@.o $(LIBRARY) -o $@$(EXT) $(FLINK_POST)


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
