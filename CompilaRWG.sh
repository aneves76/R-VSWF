#/bin/bash
gcc -std=gnu99 -I/usr/local/lib64/R/include -I/usr/local/include -fpic -g -O2 -c RectangularWaveGuide.c -o RectangularWaveGuide.o
gcc -std=gnu99 -shared -L/usr/local/lib64 -o RectangularWaveGuide.so RectangularWaveGuide.o -lm 
