#IGNORE:
CMPLR ?= gfortran
FFLAGS += -fPIC
LFLAGS = -O4

SRCE = mercury6_2.for 
DFLE = big.in  message.in  small.in element.in param.in swift.inc \
	close.in files.in mercury.inc mercury6.man  
TARFILES = $(SRCE) $(DFLE) README.txt Makefile

EXECUTABLES = mercury
OBJT = $(SRCE:.for=.o)

.for.o:
	$(CMPLR) -c $(FFLAGS) $< -DTOOLBOX

all: $(EXECUTABLES)

mercury: $(OBJT) 
	$(CMPLR) $(FFLAGS) $(OBJT) -o $@

mercury6_2.o:
	$(CMPLR) -c $(FFLAGS) $(SRCE) -c $< -DTOOLBOX

clean:
	rm -rf *.o *~ $(OUTPUTFILES) $(EXECUTABLES)

tar:
	make clean; tar cvf mercury.tar $(TARFILES); gzip mercury.tar

