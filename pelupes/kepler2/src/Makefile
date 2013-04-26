
# Makefile to build the libkepler2.a library.

# Optimization may cause failure on some systems.  Replace by -O0 in
# that case.

AMUSE_DIR?=../../../../..
OPT	  = # -O4 don't use O4 but O2 or something

GPU	  = -DGPU
CXX	  ?= mpicxx
MPICXX ?= mpicxx
CFLAGS   += $(OPT)
CXXFLAGS += $(CFLAGS) 

CODELIB = libkepler2.a

CODEOBJS = kepler-lib.o hdyn.o util.o

AR = ar ruv
RANLIB = ranlib
RM = rm

all: $(CODELIB) 

$(CODELIB): $(CODEOBJS)
	$(RM) -f $@
	$(AR) $@ $(CODEOBJS)
	$(RANLIB) $@

.cc.o: $<
	$(CXX) $(CXXFLAGS) -c -o $@ $< 

#-----------------------------------------------------------------------

kepler-lib.o:	kepler.cc
	$(MPICXX) $(CFLAGS) -c kepler.cc -o kepler-lib.o

#-----------------------------------------------------------------------

tar:
	tar czf kepler.tgz Makefile.kepler *.h *.cc

clean:
	$(RM) -r -f *.o *~ *.tgz *.a
