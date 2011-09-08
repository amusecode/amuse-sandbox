LIBNAME = libsocketrpc.a
OBJS = socketrpc.o

MPIF90 ?= openmpif90
MPICC ?= openmpicc

all: dummy_code

CFLAGS += -O2

dummy_code: dummy-code.c amuse-rpc.h $(LIBNAME)
	$(MPICC) $(CFLAGS) $(LIBNAME) -Wall -o $@ $<

$(LIBNAME): $(OBJS)
	ar -r  $(LIBNAME) $(OBJS)
	ranlib $(LIBNAME)

%.o: %.c
	$(MPICC) $(CFLAGS) -Wall -c -o $@ $<
	
%.o: %.f90
	$(MPIF90) -c -o $@ $<

clean:
	rm -Rf *.o *.lo *.a *.mo *.mod
