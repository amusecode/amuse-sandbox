all: kepler test

kepler: kepler.c
	gcc -c kepler.c -o kepler.o -lm -lmpfr -lgmp
	ar -r libkepler.a kepler.o
	ranlib libkepler.a
	mv libkepler.a ../lib

test: test.c
	gcc test.c -o test -lm -L../lib -lkepler

clean:
	rm kepler
	rm test