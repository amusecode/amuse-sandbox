
AMUSE_DIR =../../..

# fix for MacPorts+OpenMPI compilers
ifneq (,$(shell type -P openmpicxx))
	MPICXX = openmpicxx
	MPICC = openmpicc
endif

MPICXX ?= mpicxx
MPICC ?= mpicc

# -std=c99 is needed for the LOGSYS_TS macro, may be removed if deemed necessary
# -D_GNU_SOURCE is needed to define M_SQRT2 when compiling on certain machines (gpu-1@strw)
CCFLAGS += -O2 -Wall -std=c99 -D_GNU_SOURCE
CXXFLAGS += -O2 -Wall -D_GNU_SOURCE
LIBS += -lm
INCLUDE =

# set OpenCL and STDCL libraries
OPENCL = /home/inti/troep/ati-stream-sdk-v2.3-lnx64
STDCL = /home/inti/libraries/coprthr

OBJS = interface.o evolve.o 

CODE_GENERATOR = $(AMUSE_DIR)/build.py

.PHONY: makefile_debug
makefile_debug:
	@echo $(AMUSE_DIR)
	@echo $(MPICXX)
	@echo $(MPICC)

all:	huayno_worker

clean:
	rm -f *.o *.bck *.pyc *.clh worker_code.cc worker_code.h

distclean: clean
	rm -f worker worker_cl worker_mp

worker_code.cc: interface.py
	$(CODE_GENERATOR) --type=c interface.py HuaynoInterface -o $@

worker_code.h: interface.py
	$(CODE_GENERATOR) --type=h interface.py HuaynoInterface -o $@

huayno_worker: __init__.py worker_code.cc worker_code.h $(OBJS)
	$(MPICXX) $(CXXFLAGS) $(INCLUDE) worker_code.cc $(OBJS) -o $@ $(LIBS)

huayno_worker_mp: CXXFLAGS   += -fopenmp
huayno_worker_mp: __init__.py worker_code.cc worker_code.h $(OBJS)
	$(MPICXX) $(CXXFLAGS) $(INCLUDE) worker_code.cc $(OBJS) -o $@ $(LIBS)

huayno_worker_cl: CXXFLAGS   += -fopenmp -DEVOLVE_OPENCL
huayno_worker_cl: LIBS += -L$(OPENCL)/lib/x86_64 -lOpenCL -lpthread -ldl -L$(STDCL)/lib -lstdcl
huayno_worker_cl: INCLUDE += -I$(OPENCL)/include -I$(STDCL)/include
huayno_worker_cl: __init__.py worker_code.cc worker_code.h evolve_kern.clh evolve_cl.o $(OBJS) 
	$(MPICXX) $(CXXFLAGS) $(INCLUDE) worker_code.cc evolve_cl.o $(OBJS) -o $@ $(LIBS)

.cc.o: $<
	$(MPICXX) $(CXXFLAGS) $(INCLUDE) -c -o $@ $< 

.c.o: $<
	$(MPICC) $(CCFLAGS) $(INCLUDE) -c -o $@ $< 

%.clh: %.cl
	awk 'BEGIN{print "const char srcstr[]=" } {print "\""$$0"\\n\""} END{print ";"}' $< > $@
