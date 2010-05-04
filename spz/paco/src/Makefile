CSOURCE = PACO.cpp 
HEADERS = PACO.h
DATA = signal.dat README
OUTPUTFILES = ACdata.dat PACO.dat Delta.dat 

TARFILES = $(CSOURCE) $(HEADERS) $(DATA) Makefile

CXX = g++
CXXFLAGS = -Wall -O4

EXECUTABLES = paco

all: $(EXECUTABLES)

paco: PACO.o PACO.h
	$(CXX) $(CXXFLAGS) $(HEADERS) $(CSOURCE) -o $@ 

.cpp.o:
	$(CXX) -c $(CXXFLAGS) $< -DTOOLBOX

clean:
	rm -rf *.o *~ $(OUTPUTFILES) $(EXECUTABLES)

tar:
	make clean; tar cvf PACO.tar $(TARFILES); gzip PACO.tar


