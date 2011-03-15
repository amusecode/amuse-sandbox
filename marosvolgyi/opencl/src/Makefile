INCLUDE = -I"/home/marcell/downloads/ati-stream-sdk-v2.3-lnx32/ati-stream-sdk-v2.3-lnx32/include"
INCLUDEE = -I"/home/marcell/downloads/ati-stream-sdk-v2.3-lnx32/ati-stream-sdk-v2.3-lnx32/samples/opencl/SDKUtil/include/"
OPENCLLIB = -L"/home/marcell/downloads/ati-stream-sdk-v2.3-lnx32/ati-stream-sdk-v2.3-lnx32/lib/x86" -lOpenCL
OPENSDKLIB = -L"/home/marcell/downloads/ati-stream-sdk-v2.3-lnx32/ati-stream-sdk-v2.3-lnx32/samples/opencl/SDKUtil/build/debug/x86" -lSDKUtil
GLIBS = -lGL -lGLU -lglut

all: libopencl.a

libopencl.a: nbody.o
	rm -f libopencl.a
	ar -r libopencl.a nbody.o
	ranlib libopencl.a

nbody.o: NBody.cpp 
	g++ NBody.cpp $(INCLUDE) $(INCLUDEE) $(OPENCLLIB) $(OPENSDKLIB) $(GLIBS) -c -o nbody.o

clean:
	rm nbody.o
	rm libopencl.a