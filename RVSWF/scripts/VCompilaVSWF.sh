#/bin/bash
gcc -std=gnu99 -I/usr/local/lib64/R/include  -I/usr/local/include    -fpic  -g -O2 -c VVectorSphericalWaveFunctions.c -o VVectorSphericalWaveFunctions.o
gcc -std=gnu99 -shared -L/usr/local/lib64 -o VVectorSphericalWaveFunctions.so VVectorSphericalWaveFunctions.o -lm -lgsl -lgslcblas
