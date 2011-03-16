#ATISTREAMSDK = "/home/marcell/downloads/ati-stream-sdk-v2.3-lnx32/ati-stream-sdk-v2.3-lnx32"
ATISTREAMSDK = "/disks/koppoel1/marosvolgyi/troep/ati-stream-sdk-v2.3-lnx64"

INCLUDE = -I$(ATISTREAMSDK)/include
INCLUDEE = -I$(ATISTREAMSDK)/samples/opencl/SDKUtil/include
OPENCLLIB = -L$(ATISTREAMSDK)/lib/x86_64 -lOpenCL
OPENSDKLIB = -L$(ATISTREAMSDK)/samples/opencl/SDKUtil/build/debug/x86_64 -lSDKUtil
GLIBS = -lGL -lGLU -lglut

MPICXX   ?= mpicxx
CXX = $(MPICXX)

CFLAGS   += -Wall -g
CXXFLAGS += $(CFLAGS) 
LDFLAGS  += -lm $(MUSE_LD_FLAGS)

OBJS = interface.o

CODELIB = src/libopencl.a

AMUSE_DIR?=../../..

CODE_GENERATOR = $(AMUSE_DIR)/build.py

all: opencl_worker 

clean:
	$(RM) -f *.so *.o *.pyc worker_code.cc worker_code.h 
	$(RM) *~ opencl_worker worker_code.cc
	make -C src clean

$(CODELIB):
	make -C src all

worker_code.cc: interface.py
	$(CODE_GENERATOR) --type=c interface.py OpenclInterface -o $@

worker_code.h: interface.py
	$(CODE_GENERATOR) --type=H interface.py OpenclInterface -o $@

opencl_worker: worker_code.cc worker_code.h $(CODELIB) $(OBJS)
	@echo "compiling OPENCL worker"
	$(CXX) $(CXXFLAGS) $< $(OBJS) $(CODELIB) -o $@ -lm $(GLIBS) $(OPENCLLIB) $(OPENSDKLIB)
	@echo "success in compiling opencl worker"

.cc.o: $<
	@echo "compiling interface.cc"
	$(CXX) $(CXXFLAGS) -c -o $@ $<
	@echo "success in compiling interface.cc"  
