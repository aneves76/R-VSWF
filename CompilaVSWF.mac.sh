gcc -arch x86_64 -std=gnu99 -I/Library/Frameworks/R.framework/Resources/include -I/Library/Frameworks/R.framework/Resources/include/x86_64 -DNDEBUG -I/usr/local/include -fPIC  -g -O2 -c VectorSphericalWaveFunctions.c -o VectorSphericalWaveFunctions.o
gcc -arch x86_64 -std=gnu99 -dynamiclib -Wl,-headerpad_max_install_names -undefined dynamic_lookup -single_module -multiply_defined suppress -L/usr/local/lib -o VectorSphericalWaveFunctions.so VectorSphericalWaveFunctions.o -F/Library/Frameworks/R.framework/.. -framework R -Wl,-framework -Wl,CoreFoundation -lm -lgsl -lgslcblas
