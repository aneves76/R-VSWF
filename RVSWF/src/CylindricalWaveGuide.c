/* CYLINDRICAL WAVE GUIDE 
 *  *
 *  */
#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <stdlib.h>
#include <gsl/gsl_sf_bessel.h>
//------------------------------------------------------------------------------
// CYLINDRICAL WAVE GUIDE - TM MODE
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
void CWG(
      int *MD,int *M,int *S,
      double *gamma, double *kz,
      double *x,  double *y,  double *z,
      double complex *Hm, double complex *Hz, double complex *Hp,
      double complex *Em, double complex *Ez, double complex *Ep
      ){
   double k=sqrt(*gamma*(*gamma)+*kz*(*kz));
   double ctg=*kz/(*gamma);
   double csc=k/(*gamma);
   int msm=*M-*S;
   int msp=*M+*S;
   //
   *Em= PsiWG(&msm,S,gamma,kz,x,y,z)*( I*(*S)*ctg)/sqrt(2.0);
   *Ez= PsiWG(M,   S,gamma,kz,x,y,z);
   *Ep= PsiWG(&msp,S,gamma,kz,x,y,z)*(-I*(*S)*ctg)/sqrt(2.0);
   // Se TM, H <- -H, E <- E
   // Se TE, E <-  H, H <- E
   *Hm=-(*MD)*PsiWG(&msm,S,gamma,kz,x,y,z)*(*S*csc)/sqrt(2.0);
   *Hz=0;
   *Hp=-(*MD)*PsiWG(&msp,S,gamma,kz,x,y,z)*(*S*csc)/sqrt(2.0);
}
//------------------------------------------------------------------------------
// POSITION CALCULATIONS - LOOPS - COMPLETE CALCULATIONS
//------------------------------------------------------------------------------
void CylindricalWaveGuide(
                          int *TE, int *m, int *s,
                          int *nx, int *ny, int *nz,
                          double *gamma, double *kz,
                          double *x,  double *y,  double *z,
                          double *rx, double *ry, double *rz,
                          double complex *Hm, double complex *Hz, double complex *Hp,
                          double complex *Em, double complex *Ez, double complex *Ep
                    ){ 
//------------------------------------------------------------------------------
   int MD;
   int i=0;
   for(int ix=0; ix<*nx; ix++){
      for(int iy=0; iy<*ny; iy++){
         for(int iz=0; iz<*nz; iz++){
           rx[i]=x[ix];
           ry[i]=y[iy];
           rz[i]=z[iz];
           if(*TE==1){
              // MODO TE
              MD=1;
              CWG(&MD,m,s,gamma,kz,&x[ix],&y[iy],&z[iz],&Em[i],&Ez[i],&Ep[i],&Hm[i],&Hz[i],&Hp[i]);
           }
           if(*TE==0){
              // MODO TM
              MD=-1;
              CWG(&MD,m,s,gamma,kz,&x[ix],&y[iy],&z[iz],&Hm[i],&Hz[i],&Hp[i],&Em[i],&Ez[i],&Ep[i]);
           }
           i++; 
         }    
      }    
   }    
//------------------------------------------------------------------------------
}

