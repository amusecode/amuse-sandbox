ATISTREAMSDK = "/home/marcell/downloads/ati-stream-sdk-v2.3-lnx32/ati-stream-sdk-v2.3-lnx32"
#ATISTREAMSDK = "/disks/koppoel1/marosvolgyi/troep/ati-stream-sdk-v2.3-lnx64"

INCLUDE = -I$(ATISTREAMSDK)/include
INCLUDEE = -I$(ATISTREAMSDK)/samples/opencl/SDKUtil/include
OPENCLLIB = -L$(ATISTREAMSDK)/lib/x86 -lOpenCL
OPENSDKLIB = -L$(ATISTREAMSDK)/samples/opencl/SDKUtil/build/debug/x86 -lSDKUtil
GLIBS = -lGL -lGLU -lglut

all: libopencl.a

libopencl.a: nbody.o
	rm -f libopencl.a
	ar -r libopencl.a nbody.o
	ranlib libopencl.a

nbody.o: NBodylib.cpp 
	g++ NBodylib.cpp $(INCLUDE) $(INCLUDEE) $(OPENCLLIB) $(OPENSDKLIB) $(GLIBS) -c -o nbody.o

nbody: NBody.cpp
	g++ NBody.cpp $(INCLUDE) $(INCLUDEE) $(OPENCLLIB) $(OPENSDKLIB) $(GLIBS) -o nbody

clean:
	rm nbody.o
	rm libopencl.a
