CPPC=g++
CC=gcc
CFLAGS=-I/usr/share/R/include /usr/include/gsl/ -fPIC -Wall -O3 -pedantic 
LFLAGS=-L/usr/lib/gsl/ -fPIC -Wall -O3 -pedantic 

libtwovarcomp: main.o twovarcomp.o exception.o moorepenrose.o 
	${CPPC} ${LFLAGS} -shared -Wl,-soname,libtwovarcomp.so main.o twovarcomp.o exception.o moorepenrose.o -lblas -lgsl -lgslcblas -fopenmp -o libtwovarcomp.so.0.01
main.o: main.c main.h
	${CPPC} ${CFLAGS} -c main.c -o main.o -fopenmp
twovarcomp.o: twovarcomp.c twovarcomp.h exception.o moorepenrose.o
	${CPPC} ${CFLAGS} -c twovarcomp.c -o twovarcomp.o
exception.o: exception.c exception.h
	${CPPC} ${CFLAGS} -c exception.c -o exception.o
moorepenrose.o: moorepenrose.c moorepenrose.h
	${CPPC} ${CFLAGS} -c moorepenrose.c -o moorepenrose.o
clean: 
	rm -f *.o *.so *~ *.so.*
