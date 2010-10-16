all: kepler test demo

kepler: kepler.c
	gcc -c kepler.c -o kepler.o -lm -lmpfr -lgmp
	ar -r libkepler.a kepler.o
	ranlib libkepler.a
	mv libkepler.a ../lib

test: test.c
	gcc test.c -o test -lm -lmpfr -lgmp -L../lib -lkepler 

demo: demo.c
	gcc demo.c -o demo -lm -lmpfr -lgmp -L../lib -lkepler -lSDL -L../lib -ldraw

clean:
	rm kepler.o
	rm ../lib/libkepler.a
	rm test
	rm demo