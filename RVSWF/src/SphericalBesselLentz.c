/*------------------------------------------------------------------------------*/
#include <stdio.h>
#include <math.h>
/*------------------------------------------------------------------------------*/
double Tk(double x, int k){
   return((2*k+1)/x);
}
/*------------------------------------------------------------------------------*/
/* Spherical Bessel Functions calculation via Lentz method                      */
void SphericalBesselLentz(int *n, double *x, double *jl, double *djl){
   if(*n<2*(int *xo)){
      n=2*(int *xo);
   }
   n=*n+10;
}
