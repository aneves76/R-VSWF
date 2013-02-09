#/bin/bash
gcc -std=gnu99 -I/usr/local/lib64/R/include -I/usr/local/include -fpic -g -O2 -c CylindricalWaveGuide.c -o CylindricalWaveGuide.o
gcc -std=gnu99 -shared -L/usr/local/lib64 -o CylindricalWaveGuide.so CylindricalWaveGuide.o -lm -lgsl -lgslcblas
