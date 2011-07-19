CFLAGS   += -Wall -g
CXXFLAGS += $(CFLAGS) 
LDFLAGS  += -lm $(MUSE_LD_FLAGS)

CODELIB = libclustorvms.a
OBJS = integrator_MS
CODEOBJS = test.o

AR = ar ruv
RANLIB = ranlib
RM = rm

all: $(OBJS) $(CODELIB) 

clean:
	$(RM) -f *.o *.a

$(CODELIB): $(CODEOBJS)
	$(RM) -f $@
	$(AR) $@ $(CODEOBJS)
	$(RANLIB) $@

.cc.o: $<
	$(CXX) $(CXXFLAGS) -c -o $@ $< 

# path to libmpfr.a and libgmp.a (platform dependant)
LDIR = ./Packages/
# path to dlmalloc.o, dlmalloc.c, mpreal.cpp, mpreal.h
LLDIR = ./Packages/mpfr_c++/

###################################################################

C   = g++ -O2 -c
CC  = g++ -O2 -o

CP  = g++ -O2 -pg -c
CCP = g++ -O2 -pg -o

###################################################################

DIR_int_ms = ./

###################################################################

###################################################################
## Run command
###################################################################

runMS:
	 $(DIR_int_ms)./integrator_MS.exe

###################################################################
## Compile commands
###################################################################

integrator_MS: $(DIR_int_ms)STAR.h $(DIR_int_ms)FORCE.h $(DIR_int_ms)DYNAMICS.h $(DIR_int_ms)CLOCK.h $(DIR_int_ms)BS_INTEGRATOR.h $(DIR_int_ms)CLUSTER.h mpreal.o dlmalloc.o
	$(CC) $(DIR_int_ms)integrator_MS.exe $(DIR_int_ms)integrator_MS.cpp $(DIR_int_ms)mpreal.o $(DIR_int_ms)dlmalloc.o $(LDIR)libmpfr.a $(LDIR)libgmp.a

###################################################################

mpreal.o: $(LLDIR)mpreal.cpp $(LLDIR)mpreal.h
	$(C) $(LLDIR)mpreal.cpp
	#mv mpreal.o $(DIR_int_ms)

dlmalloc.o: $(LLDIR)dlmalloc.c $(LLDIR)dlmalloc.h
	$(C) $(LLDIR)dlmalloc.c
	#mv dlmalloc.o $(DIR_int_ms)

test.o: 
	mpicc -Wall -g -c test.cc -o test.o -lm $(LDIR)libmpfr.a $(LDIR)libgmp.a

###################################################################

clean:
	rm -f  $(DIR_int_ms)$(OBJS) *~ $(DIR_int_ms)*~  $(DIR_int_ms)*.o $(DIR_int_ms)*.exe


