all: libdraw.a

libdraw.o: draw.c 
	gcc -c draw.c -o libdraw.o -lSDL -lm

libdraw.a: libdraw.o
	ar -r libdraw.a libdraw.o
	ranlib libdraw.a
	mv libdraw.a ../lib

clean: 
	rm *.o
	rm *.a