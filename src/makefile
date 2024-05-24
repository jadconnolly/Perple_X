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
#                  COMP77 = the name of the local fortran 77 compatible compiler

COMP77 = gfortran

#                  Architecture

FARCH = -m64

#                  FFLAGS = compiler options
#                  FLINK = linker options
# -Wpedantic 
# -Wunused
# -finit-local-zero 
# -ffpe-trap=zero,overflow,underflow to catch fp errors.
# -ffpe-trap 
# -fbounds-check array bound checking
# -finit-real=snan -fpe=trap check for unitialized variables
# -Wall -Wextra for extensive warnings
# -static-libgfortran
# -m64, -m32  
# -Os, -Ofast optimize size, speed

# For distribution using static (pre-linked) compiler run-time libraries, use
# following compiler flags for gfortran builds on MacOSX/Darwin.  gfortran
# versions up through 4.9 lack the ability to include a static libquadmath, so
# awkward link flags (FLINK) and directory for static libquadmath defined.
# Most commonly, gfortran will use the directory /usr/local/lib, but some
# builds will differ on location.  To find your library location, use a command
# (from the command line) like,
#    find /usr/local -name libquadmath.a -print
# and substitute the directory name you find for LQMDIR.  En garde choosing
# 32/64 bit libraries.
# G. Helffrich/30 Nov. '15

# George's justifications/objections for/to current f77 options:

# -fPIC (position independent code) added for perplexwrap, there may be some OS / machine 
# combinations where -fPIC won’t work, but probably the Julia linkage won’t work on those 
# either, so no functionality will be lost.  The Julia wrappers shouldn’t really use -fPIC 
# explicitly; when Julia tries to link with the wrappers, then if PIC is not available, 
# they won’t link and it will be Julia’s problem to solve - not perplex's.

# -std=legacy prevents gfortran from warning about things like not ending numbered DO loops 
# with CONTINUE statements, labels on IF statements, 3-way IF statements, etc.  There are a 
# lot of these warnings that the perplex programs elicit.

# -fno-backtrace When you interrupt or a program ends with a runtime error, the Fortran 
# library gives a call backtrace.  -fno-backtrace suppresses this. The backtrace is useless 
# because it only contains hex addresses, and the compiler never gives you any info about 
# where the program was loaded anyway, or what addresses correspond to what subprograms, so it 
# is better left out. (True, unless you compile with -g, in which case the backtrace is useful) 

#FFLAGS = -O3 #original github make
#FFLAGS = -O3 -std=legacy -fno-backtrace $(FARCH) #for perplex
FFLAGS = -O3 -fPIC -Wpedantic -ffixed-line-length-none $(FARCH) #for perplex + perplexwrap

FLINK = -static-libgfortran -lgfortran -lgcc -lm $(FARCH)
LLINK = -shared -lgfortran -lgcc -lm $(FARCH)

MYOBJ = actcor build fluids ctransf frendly meemum convex pstable pspts psvdraw pssect pt2curv vertex werami MC_fit perplexwrap
all: $(MYOBJ)

LIBRARY = pscom.o pslib.o rlib.o tlib.o flib.o olib.o resub.o minime_blas.o blas2lib.o

LIBRARY_PS = cont_lib.o pscom.o pslib.o tlib.o

clean: 
	rm -f *.o $(MYOBJ) *.so *.mod

######### PerplexWrap shared library compilation #########################
# PerplexWrap won't compile under -std=legacy, and needs link option -shared?

perplexwrap: perplexwrap.o $(LIBRARY)
	$(COMP77) $(FFLAGS) $(LLINK) perplexwrap.o $(LIBRARY) -o perplexwrap.so
perplexwrap.o: perplexwrap.f
	$(COMP77) $(FFLAGS) -c perplexwrap.f

###################### TARGETS FOR 692+ PROGRAMS #########################

MC_fit: MC_fit.o $(LIBRARY)
	$(COMP77) $(FFLAGS) $(FLINK) MC_fit.o $(LIBRARY) -o $@$(EXT)

convex: convex.o $(LIBRARY)
	$(COMP77) $(FFLAGS) $(FLINK) convex.o $(LIBRARY) -o $@$(EXT)

meemum: meemum.o $(LIBRARY)
	$(COMP77) $(FFLAGS) $(FLINK) meemum.o $(LIBRARY) -o $@$(EXT)

vertex: vertex.o $(LIBRARY)
	$(COMP77) $(FFLAGS) $(FLINK) vertex.o $(LIBRARY) -o $@$(EXT)

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

#hptover: hptover.o
#	$(COMP77) $(FFLAGS) $@.o -o $@
#htog: htog.o
#	$(COMP77) $(FFLAGS) $@.o -o $@
# targets missing from '07:
#rk: rk.o flib.o tlib.o
#	$(COMP77) $(FFLAGS) $@.o tlib.o flib.o -o $@
#ps_p_contor: ps_p_contor.o pslib.o
#	$(COMP77) $(FFLAGS) $@.o pslib.o -o $@
#ge0pt: ge0pt.o
#	$(COMP77) $(FFLAGS) $@.o -o $@
#satsurf: satsurf.o flib.o tlib.o
#	$(COMP77) $(FFLAGS) $@.o flib.o tlib.o -o $@
#rewrite,sox,gox,cohscont 

#################################################################################
# Explicitly make objects to override the default compiler choice implemented
# on some machines.

actcor.o: actcor.f perplex_parameters.h
	$(COMP77) $(FFLAGS) -c actcor.f
build.o: build.f perplex_parameters.h
	$(COMP77) $(FFLAGS) -c build.f
DEW_2_ver.o: DEW_2_ver.f perplex_parameters.h
	$(COMP77) $(FFLAGS) -c DEW_2_ver.f
fluids.o: fluids.f perplex_parameters.h
	$(COMP77) $(FFLAGS) -c fluids.f
convex.o: convex.f perplex_parameters.h
	$(COMP77) $(FFLAGS) -c convex.f
cont_lib.o: cont_lib.f perplex_parameters.h
	$(COMP77) $(FFLAGS) -c cont_lib.f
ctransf.o: ctransf.f perplex_parameters.h
	$(COMP77) $(FFLAGS) -c ctransf.f
frendly.o: frendly.f perplex_parameters.h
	$(COMP77) $(FFLAGS) -c frendly.f
#ge0pt.o: ge0pt.f perplex_parameters.h
#	$(COMP77) $(FFLAGS) -c ge0pt.f
#hptover.o: hptover.f perplex_parameters.h
#	$(COMP77) $(FFLAGS) -c hptover.f
#htog.o: htog.f perplex_parameters.h
#	$(COMP77) $(FFLAGS) -c htog.f
psect.o: psect.f perplex_parameters.h
	$(COMP77) $(FFLAGS) -c psect.f
psvdraw.o: psvdraw.f perplex_parameters.h
	$(COMP77) $(FFLAGS) -c psvdraw.f
pscom.o: pscom.f perplex_parameters.h
	$(COMP77) $(FFLAGS) -c pscom.f
pspts.o: pspts.f perplex_parameters.h
	$(COMP77) $(FFLAGS) -c pspts.f
pstable.o: pstable.f perplex_parameters.h
	$(COMP77) $(FFLAGS) -c pstable.f
pt2curv.o: pt2curv.f perplex_parameters.h
	$(COMP77) $(FFLAGS) -c pt2curv.f
werami.o: werami.f perplex_parameters.h
	$(COMP77) $(FFLAGS) -c werami.f
flib.o: flib.f perplex_parameters.h
	$(COMP77) $(FFLAGS) -c flib.f
pslib.o: pslib.f perplex_parameters.h
	$(COMP77) $(FFLAGS) -c pslib.f
vertex.o: vertex.f perplex_parameters.h
	$(COMP77) $(FFLAGS) -c vertex.f
meemum.o: meemum.f perplex_parameters.h
	$(COMP77) $(FFLAGS) -c meemum.f
MC_fit.o: MC_fit.f perplex_parameters.h
	$(COMP77) $(FFLAGS) -c MC_fit.f
resub.o: resub.f perplex_parameters.h
	$(COMP77) $(FFLAGS) -c resub.f
tlib.o: tlib.f perplex_parameters.h
	$(COMP77) $(FFLAGS) -c tlib.f
olib.o: olib.f perplex_parameters.h
	$(COMP77) $(FFLAGS) -c olib.f
rlib.o: rlib.f perplex_parameters.h
	$(COMP77) $(FFLAGS) -c rlib.f
blas2lib.o: blas2lib.f perplex_parameters.h
	$(COMP77) $(FFLAGS) -c blas2lib.f
minime_blas.o: minime_blas.f perplex_parameters.h
	$(COMP77) $(FFLAGS) -c minime_blas.f
