/* CYLINDRICAL WAVE GUIDE 
 *  *
 *  */
#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <stdlib.h>
#include <gsl/gsl_sf_bessel.h>
//------------------------------------------------------------------------------
double complex PsiWG(int *m,int *s,double *gamma,double *kz,double *x,double *y,double *z){
   double rho=sqrt(*x*(*x)+*y*(*y));
   double cph,sph;
   double zero=1e-15;
   cph=*x/rho;
   sph=*y/rho;
   double complex eiph=cpow(cph+I*sph,*m*(*s));
   double complex eikz=cexp(I*(*kz)*(*z)); 
   double jn=gsl_sf_bessel_Jn(*m,rho*(*gamma));
   double complex u=jn*eiph*eikz;
   return(u);
}
//------------------------------------------------------------------------------
void BBP(
      int *P,int *M,int *S,
      double *gamma, double *kz,
      double *x,  double *y,  double *z,
      double complex *Hm, double complex *Hz, double complex *Hp,
      double complex *Em, double complex *Ez, double complex *Ep
      ){
   double k=sqrt(*gamma*(*gamma)+*kz*(*kz));
   double cth=*kz/k;
   double sth=(*gamma)/k;

   int pp1=1+*P;
   int pm1=1-*P;
   int mm1=*M-1;
   int mspm1=*M-1+*S*(*P-1);
   int mspm0=*M-1+*S*(*P);
   int mspp1=*M-1+*S*(*P+1);

   *Em=pp1*PsiWG(&mm1,S,gamma,kz,x,y,z)/2-*P*(sth*sth)*PsiWG(&mspm1,S,gamma,kz,x,y,z)/2;
   *Ez= -I*(*P)*(*S)*cth*sth*PsiWG(&mspm0,S,gamma,kz,x,y,z)/sqrt(2.0);
   *Ep=pm1*PsiWG(&mm1,S,gamma,kz,x,y,z)/2+*P*(sth*sth)*PsiWG(&mspp1,S,gamma,kz,x,y,z)/2;
   *Hm=-I*(*P)*cth*PsiWG(&mm1,S,gamma,kz,x,y,z)*pp1/2.0;
   *Hz=-(*S)*sth*PsiWG(&mspm0,S,gamma,kz,x,y,z)/sqrt(2.0);
   *Hp=-I*(*P)*cth*PsiWG(&mm1,S,gamma,kz,x,y,z)*pm1/2.0;
}
//------------------------------------------------------------------------------
// POSITION CALCULATIONS - LOOPS - COMPLETE CALCULATIONS
//------------------------------------------------------------------------------
void BesselBeamsP(
                          int *p, int *m, int *s, 
                          int *nx, int *ny, int *nz,
                          double *gamma, double *kz,
                          double *x,  double *y,  double *z,
                          double *rx, double *ry, double *rz,
                          double complex *Hm, double complex *Hz, double complex *Hp,
                          double complex *Em, double complex *Ez, double complex *Ep
                    ){ 
//------------------------------------------------------------------------------
   int i=0;
   for(int ix=0; ix<*nx; ix++){
      for(int iy=0; iy<*ny; iy++){
         for(int iz=0; iz<*nz; iz++){
           rx[i]=x[ix];
           ry[i]=y[iy];
           rz[i]=z[iz];
           BBP(p,m,s,gamma,kz,&x[ix],&y[iy],&z[iz],&Hm[i],&Hz[i],&Hp[i],&Em[i],&Ez[i],&Ep[i]);
           i++; 
         }    
      }    
   }    
//------------------------------------------------------------------------------
}

