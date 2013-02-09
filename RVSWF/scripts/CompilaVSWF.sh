#/bin/bash
gcc -std=gnu99 -I/usr/local/lib64/R/include  -I/usr/local/include    -fpic  -g -O2 -c VectorSphericalWaveFunctions.c -o VectorSphericalWaveFunctions.o
gcc -std=gnu99 -shared -L/usr/local/lib64 -o VectorSphericalWaveFunctions.so VectorSphericalWaveFunctions.o -lm -lgsl -lgslcblas
