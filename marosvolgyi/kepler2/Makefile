all: kepler test

kepler: kepler.c
	gcc -c kepler.c -o kepler.o -lm -lmpfr -lgmp
	ar -r libkepler.a kepler.o
	ranlib libkepler.a

test: test.c
	gcc test.c -o test -lm -L./ -lkepler

clean:
	rm -f *.o *.a *.dat test
