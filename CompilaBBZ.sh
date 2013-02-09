#/bin/bash
gcc -std=gnu99 -I/usr/local/lib64/R/include -I/usr/local/include -fpic -g -O2 -c BesselBeamsZ.c -o BesselBeamsZ.o
gcc -std=gnu99 -shared -L/usr/local/lib64 -o BesselBeamsZ.so BesselBeamsZ.o -lm -lgsl -lgslcblas
