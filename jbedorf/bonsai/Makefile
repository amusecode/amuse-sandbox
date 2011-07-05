MPICXX   ?= mpicxx

CFLAGS   += -Wall -g
CXXFLAGS += $(CFLAGS) 
LDFLAGS  += -lm $(MUSE_LD_FLAGS)

CUDA_TK = /usr/local/cuda/

CXXFLAGS += -I$(CUDA_TK)/include
LIBS = -lcuda -lOpenCL

OBJS = interface.o libraryInterface.o

CODELIB = src/libbonsai.a

CODELIB += $(LIBS)

AMUSE_DIR?=../../../..

CODE_GENERATOR = $(AMUSE_DIR)/build.py

all: bonsai_worker 

clean:
	$(RM) -f *.so *.o *.pyc worker_code.cc worker_code.h 
	$(RM) *~ bonsai_worker worker_code.cc
	make -C src clean

$(CODELIB):
	make -C src all

worker_code.cc: interface.py
	$(CODE_GENERATOR) --type=c interface.py BonsaiInterface -o $@

worker_code.h: interface.py
	$(CODE_GENERATOR) --type=H interface.py BonsaiInterface -o $@

bonsai_worker: worker_code.cc worker_code.h $(CODELIB) $(OBJS)
	$(MPICXX) $(CXXFLAGS) $< $(OBJS) $(CODELIB) -o $@

.cpp.o: $<
	$(MPICXX) $(CXXFLAGS) -c -o $@ $< 
.cc.o: $<
	$(MPICXX) $(CXXFLAGS) -c -o $@ $< 
