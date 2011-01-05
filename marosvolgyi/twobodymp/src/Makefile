all: kepler

kepler: kepler.c kepler.h
	mpicc -Wall -g -c kepler.c -o kepler.o -lm -lmpfr -lgmp
	ar -r libkepler.a kepler.o
	ranlib libkepler.a

demo: demo.c
	gcc demo.c -o demo -lm -lmpfr -lgmp -L./ -lkepler 

clean:
	rm kepler.o
	rm libkepler.a
