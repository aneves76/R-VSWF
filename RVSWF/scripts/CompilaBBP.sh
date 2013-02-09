#/bin/bash
gcc -std=gnu99 -I/usr/local/lib64/R/include -I/usr/local/include -fpic -g -O2 -c BesselBeamsP.c -o BesselBeamsP.o
gcc -std=gnu99 -shared -L/usr/local/lib64 -o BesselBeamsP.so BesselBeamsP.o -lm -lgsl -lgslcblas
