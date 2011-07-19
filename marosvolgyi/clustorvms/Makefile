MPICXX   ?= mpicxx

CFLAGS   += -Wall -g
CXXFLAGS += $(CFLAGS) 
LDFLAGS  += -lm $(MUSE_LD_FLAGS)

OBJS = interface.o

CODELIB = ./src/libclustorvms.a

AMUSE_DIR?=../../../

CODE_GENERATOR = $(AMUSE_DIR)/build.py

# path to libmpfr.a and libgmp.a (platform dependant)
LDIR = ./src/Packages/
# path to dlmalloc.o, dlmalloc.c, mpreal.cpp, mpreal.h
LLDIR = ./src/Packages/mpfr_c++/
DIR_int_ms = ./src/

all: tests clustorvms_worker

tests:  
	make -C src all

clean:
	$(RM) -f *.so *.o *.pyc worker_code.cc worker_code.h 
	$(RM) *~ clustorvms_worker worker_code.cc
	make -C src clean

$(CODELIB):
	make -C src all

worker_code.cc: interface.py
	$(CODE_GENERATOR) --type=c interface.py ClustorvmsInterface -o $@

worker_code.h: interface.py
	$(CODE_GENERATOR) --type=H interface.py ClustorvmsInterface -o $@

clustorvms_worker: worker_code.cc worker_code.h $(CODELIB) $(OBJS)
	$(MPICXX) $(CXXFLAGS) $< $(OBJS) -o $@ $(CODELIB) ./src/mpreal.o $(DIR_int_ms)dlmalloc.o $(LDIR)libmpfr.a $(LDIR)libgmp.a

.cc.o: $<
	$(CXX) $(CXXFLAGS) -c -o $@ $< 
