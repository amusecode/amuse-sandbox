MPICXX   ?= mpicxx
CXX = $(MPICXX)

CFLAGS   += -Wall -g
CXXFLAGS += $(CFLAGS) 
LDFLAGS  += -lm $(MUSE_LD_FLAGS)

OBJS = interface.o

CODELIB = src/libkepler.a

AMUSE_DIR?=../../..

CODE_GENERATOR = $(AMUSE_DIR)/build.py

all: twobodymp_worker 

clean:
	$(RM) -f *.so *.o *.pyc worker_code.cc worker_code.h 
	$(RM) *~ twobodymp_worker worker_code.cc
	make -C src clean

$(CODELIB):
	make -C src all

worker_code.cc: interface.py
	$(CODE_GENERATOR) --type=c interface.py twobodympInterface -o $@

worker_code.h: interface.py
	$(CODE_GENERATOR) --type=H interface.py twobodympInterface -o $@

twobodymp_worker: worker_code.cc worker_code.h $(CODELIB) $(OBJS)
	@echo "compiling twobodymp worker"
	$(CXX) $(CXXFLAGS) $< $(OBJS) $(CODELIB) -o $@ -lm -lmpfr -lgmp -lSDL -L../lib -ldraw
	@echo "success in compiling twobodymp worker"

.cc.o: $<
	@echo "compiling interface.cc"
	$(CXX) $(CXXFLAGS) -c -o $@ $<
	@echo "success in compiling interface.cc"  
