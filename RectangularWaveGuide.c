/* RECTANGULAR WAVE GUIDE */
#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <stdlib.h>
//------------------------------------------------------------------------------
// RECTANGULAR WAVE GUIDE - TM MODE
//------------------------------------------------------------------------------
void RWG_TM(
      double *kx, double *ky, double *kz,
      double *x,  double *y,  double *z,
      double complex *Hm, double complex *Hz, double complex *Hp,
      double complex *Em, double complex *Ez, double complex *Ep
      ){
   double gamma2=*kx*(*kx)+*ky*(*ky);
   double k=sqrt(*kx*(*kx)+*ky*(*ky)+*kz*(*kz));
   //
   double complex ETMx=I*(*kz*(*kx)/gamma2)*cos(*kx*(*x))*sin(*ky*(*y))*cexp(I*(*kz)*(*z));
   double complex ETMy=I*(*kz*(*ky)/gamma2)*sin(*kx*(*x))*cos(*ky*(*y))*cexp(I*(*kz)*(*z));
   //
   double complex ETMm=(ETMx-I*ETMy)/sqrt(2.0);
   double complex ETMz=sin(*kx*(*x))*sin(*ky*(*y))*cexp(I*(*kz)*(*z));
   double complex ETMp=(ETMx+I*ETMy)/sqrt(2.0);
   //
   double complex HTMx=-I*(k*(*ky)/gamma2)*sin(*kx*(*x))*cos(*ky*(*y))*cexp(I*(*kz)*(*z));
   double complex HTMy= I*(k*(*kx)/gamma2)*cos(*kx*(*x))*sin(*ky*(*y))*cexp(I*(*kz)*(*z));
   //
   double complex HTMm=(HTMx-I*HTMy)/sqrt(2.0);
   double complex HTMz=0.0+I*0.0;
   double complex HTMp=(HTMx+I*HTMy)/sqrt(2.0);
   //
   Em=&ETMm;
   Ez=&ETMz;
   Ep=&ETMp;
   //
   Hm=&HTMm;
   Hz=&HTMz;
   Hp=&HTMp;
}
//------------------------------------------------------------------------------
// RECTANGULAR WAVE GUIDE - TE MODE
//------------------------------------------------------------------------------
void RWG_TE(
      double *kx, double *ky, double *kz,
      double *x,  double *y,  double *z,
      double complex *Hm, double complex *Hz, double complex *Hp,
      double complex *Em, double complex *Ez, double complex *Ep
      ){
   double gamma2=*kx*(*kx)+*ky*(*ky);
   double k=sqrt(*kx*(*kx)+*ky*(*ky)+*kz*(*kz));
   //
   double complex ETEx=-I*(k*(*ky)/gamma2)*cos(*kx*(*x))*sin(*ky*(*y))*cexp(I*(*kz)*(*z));
   double complex ETEy= I*(k*(*kx)/gamma2)*sin(*kx*(*x))*cos(*ky*(*y))*cexp(I*(*kz)*(*z));
   //
   *Em=(ETEx-I*ETEy)/sqrt(2.0);
   *Ez=0.0+I*0.0;
   *Ep=(ETEx+I*ETEy)/sqrt(2.0);
   //
   double complex HTEx=-I*(*kz*(*kx)/gamma2)*sin(*kx*(*x))*cos(*ky*(*y))*cexp(I*(*kz)*(*z));
   double complex HTEy=-I*(*kz*(*ky)/gamma2)*cos(*kx*(*x))*sin(*ky*(*y))*cexp(I*(*kz)*(*z));
   //
   *Hm=(HTEx-I*HTEy)/sqrt(2.0);
   *Hz=cos(*kx*(*x))*cos(*ky*(*y))*cexp(I*(*kz)*(*z));
   *Hp=(HTEx+I*HTEy)/sqrt(2.0);
}
//------------------------------------------------------------------------------
// POSITION CALCULATIONS - LOOPS - COMPLETE CALCULATIONS
//------------------------------------------------------------------------------
void RectangularWaveGuide(
                          int *TE,
                          int *nx, int *ny, int *nz,
                          double *kx, double *ky, double *kz,
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
           if(*TE==1){
              RWG_TE(kx,ky,kz,&x[ix],&y[iy],&z[iz],&Hm[i],&Hz[i],&Hp[i],&Em[i],&Ez[i],&Ep[i]);
           }
           if(*TE==0){
              RWG_TM(kx,ky,kz,&x[ix],&y[iy],&z[iz],&Hm[i],&Hz[i],&Hp[i],&Em[i],&Ez[i],&Ep[i]);
           }
           i++; 
         }    
      }    
   }    
//------------------------------------------------------------------------------
}

