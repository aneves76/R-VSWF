/* Mm.c
 *  *
 *  * Compile
 *  gcc -std=gnu99 -I/usr/local/lib64/R/include  -I/usr/local/include    -fpic  -g -O2 -c PHansenMultipoles.c -o PHansenMultipoles.o
 *  gcc -std=gnu99 -shared -L/usr/local/lib64 -o PHansenMultipoles.so PHansenMultipoles.o -lm -lgsl -lgslcblas
 *
 *  * gcc -fPIC -O2 -g -std=gnu99 -m64 -c Mm.c -o Mm.o
 *  * gcc -shared -Wl,-O1 -o Mm.so Mm.o -lR -lm -lgsl -lgslcblas
 *  *
 *  */
#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <stdlib.h>
#include <gsl/gsl_sf_bessel.h>
//------------------------------------------------------------------------------
// CONSTANTES PARA Qlm
//------------------------------------------------------------------------------
double alfaQ(int lo, int mo){
   double l=lo;
   double m=mo;
   return sqrt(((2*l-1)*(2*l+1))/((l-m)*(l+m)));
}
//------------------------------------------------------------------------------
double betaQ(int lo, int mo){
   double l=lo;
   double m=mo;
   return sqrt((2*l+1)/(2*l-3))*sqrt(((l+m-1)*(l-m-1))/((l-m)*(l+m)));
}
//------------------------------------------------------------------------------
double gammaQ(int lo){
   double l=lo;
   return sqrt((2*l+1)/(2*l));
}
//------------------------------------------------------------------------------
double deltaQ(int lo){
   double l=lo;
   return sqrt(2*l+1);
}
//------------------------------------------------------------------------------
// POSICIONAMENTO GERAL DOS ELEMENTOS EM RELACAO A (l,m) DOS MULTIPOLOS
// A lista em C comeca com o elemento 0
//------------------------------------------------------------------------------
int jlm(int l, int m){
   if(abs(m)>l){
      return 0;
   }else{
      return l*(l+1)+m;
   }
}
//------------------------------------------------------------------------------
// CONSTANTES PARA X
//------------------------------------------------------------------------------
double cp(int lo, int mo){
   double l=lo;
   double m=mo;
   return sqrt(l*(l+1)-m*(m+1));
}
//------------------------------------------------------------------------------
double cm(int lo, int mo){
   double l=lo;
   double m=mo;
   return sqrt(l*(l+1)-m*(m-1));
}
//------------------------------------------------------------------------------
// POSITION DEPENDENT CALCULATIONS - MULTIPLES LOOPS - COMPLETE CALCULATIONS
//------------------------------------------------------------------------------
void VectorSphericalWaveFunctionsKernel(double *k,double *x, double *y, double *z,
                    int *lmax, 
                    double complex *GTE, double complex *GTM,
                    double complex *Em, double complex *Ez, double complex *Ep,
                    double complex *Hm, double complex *Hz, double complex *Hp
                    ){
//------------------------------------------------------------------------------
   int l,m;
   int LMAX=*lmax*(*lmax+2);
   int LMAXE=(*lmax+1)*(*lmax+3);
   double cph; 
   double sph;
   double rho=sqrt(*x*(*x)+*y*(*y));
   double r=sqrt(rho*rho+*z*(*z));
   double sth=rho/r;
   double cth=*z/r;
//------------------------------------------------------------------------------
   if((*x==0)&&(*y==0)){
      cph=1;
      sph=0;
   }else{
      cph=*x/rho;
      sph=*y/rho;
   }
//------------------------------------------------------------------------------
   // Spherical Bessel Funtions
   double JLM[*lmax+2];
   gsl_sf_bessel_jl_steed_array(*lmax+1,*k*r,JLM);
//------------------------------------------------------------------------------
   // Qlm - primeiros 4 termos
   double Qlm[LMAXE];
   Qlm[jlm(0, 0)]=1/sqrt(4*M_PI);
   Qlm[jlm(1, 1)]=-gammaQ(1)*sth*Qlm[jlm(0,0)]; // Q11
   Qlm[jlm(1, 0)]=sqrt(3.0)*cth*Qlm[jlm(0,0)];  // Q10
   Qlm[jlm(1,-1)]=-Qlm[jlm(1,1)];               // Q11*(-1)
//------------------------------------------------------------------------------
   // Complex Exponencial for m=-1,0,1
   double complex Eim[2*(*lmax)+3];
   Eim[*lmax-1]=(cph-I*sph);
   Eim[*lmax  ]=1+I*0;
   Eim[*lmax+1]=(cph+I*sph);
//------------------------------------------------------------------------------
   // Ylm - primeiros 4 termos
   double complex Ylm[LMAXE];
   Ylm[jlm(0, 0)]=Qlm[jlm(0, 0)];
   Ylm[jlm(1,-1)]=Qlm[jlm(1,-1)]*Eim[*lmax-1];
   Ylm[jlm(1, 0)]=Qlm[jlm(1, 0)];
   Ylm[jlm(1, 1)]=Qlm[jlm(1, 1)]*Eim[*lmax+1];
//------------------------------------------------------------------------------
   // VECTOR SPHERICAL HARMONICS
//------------------------------------------------------------------------------
   // r
   double complex rm=sth*(cph-I*sph)/sqrt(2);
   double complex rz=cth;
   double complex rp=sth*(cph+I*sph)/sqrt(2);
   // X
   double complex XM;
   double complex XZ;
   double complex XP;
   // Y
   double complex YM;
   double complex YZ;
   double complex YP;
   // V
   double complex VM;
   double complex VZ;
   double complex VP;
//------------------------------------------------------------------------------
   // HANSEN MULTIPOLES
//------------------------------------------------------------------------------
   // M
   double complex MM;
   double complex MZ;
   double complex MP;
   // N
   double complex NM;
   double complex NZ;
   double complex NP;
//------------------------------------------------------------------------------
   // OTHERS
//------------------------------------------------------------------------------
   double kl;
//------------------------------------------------------------------------------
   // MAIN LOOP
   for(l=1;l<=(*lmax);l++){
//------------------------------------------------------------------------------
      //Qlm extremos positivos um passo a frente
      Qlm[jlm(l+1, l+1)]=-gammaQ(l+1)*sth*Qlm[jlm(l,l)];
      Qlm[jlm(l+1, l  )]= deltaQ(l+1)*cth*Qlm[jlm(l,l)];
      //Qlm extremos negativos um passo a frente
      Qlm[jlm(l+1,-l-1)]=pow(-1,l+1)*Qlm[jlm(l+1, l+1)];
      Qlm[jlm(l+1,-l  )]=pow(-1,l  )*Qlm[jlm(l+1, l  )];
      // Exponenciais um passo a frente
      Eim[*lmax+l+1]=Eim[*lmax+l]*(cph+I*sph);
      Eim[*lmax-l-1]=Eim[*lmax-l]*(cph-I*sph);
      // Harmonicos esfericos extremos um passo a frente
      Ylm[jlm(l+1, l+1)]=Qlm[jlm(l+1, l+1)]*Eim[*lmax+l+1];
      Ylm[jlm(l+1, l  )]=Qlm[jlm(l+1, l  )]*Eim[*lmax+l  ];
      Ylm[jlm(l+1,-l-1)]=Qlm[jlm(l+1,-l-1)]*Eim[*lmax-l-1];
      Ylm[jlm(l+1,-l  )]=Qlm[jlm(l+1,-l  )]*Eim[*lmax-l  ];
      // others
      kl=1/(sqrt(l*(l+1)));
//------------------------------------------------------------------------------
      for(int m=l; m>=(-l); m--){
         // CALCULATIONS OF SSH
         if(m>=0){
            Qlm[jlm(l+1, m)]=alfaQ(l+1,m)*cth*Qlm[jlm(l,m)]-betaQ(l+1,m)*Qlm[jlm(l-1,m)];
            Qlm[jlm(l+1,-m)]=pow(-1,m)*Qlm[jlm(l+1, m)];
            Ylm[jlm(l+1, m)]=Qlm[jlm(l+1, m)]*Eim[*lmax+m];
            Ylm[jlm(l+1,-m)]=Qlm[jlm(l+1,-m)]*Eim[*lmax-m];
         }
         // CALCULATIONS OF VSH
         // X
         XM=kl*cm(l,m)*Ylm[jlm(l,m-1)]/sqrt(2);
         XZ=kl*m*Ylm[jlm(l,m  )];
         XP=kl*cp(l,m)*Ylm[jlm(l,m+1)]/sqrt(2);
         // Y
         YM=rm*Ylm[jlm(l,m)];
         YZ=rz*Ylm[jlm(l,m)];
         YP=rp*Ylm[jlm(l,m)];
         // V
         VM=rm*XZ-rz*XM;
         VZ=rp*XM-rm*XP;
         VP=rz*XP-rp*XZ;
         // CALCULATION OF HANSEM MULTIPOLES
         // M
         MM=JLM[l]*XM;
         MZ=JLM[l]*XZ;
         MP=JLM[l]*XP;
         // N
         NM=((1.*l+1.)*JLM[l-1]-l*JLM[l+1])*VM/(2.*l+1.)+sqrt(l*(l+1.))*(JLM[l-1]+JLM[l+1])*YM/(2.*l+1.);
         NZ=((1.*l+1.)*JLM[l-1]-l*JLM[l+1])*VZ/(2.*l+1.)+sqrt(l*(l+1.))*(JLM[l-1]+JLM[l+1])*YZ/(2.*l+1.);
         NP=((1.*l+1.)*JLM[l-1]-l*JLM[l+1])*VP/(2.*l+1.)+sqrt(l*(l+1.))*(JLM[l-1]+JLM[l+1])*YP/(2.*l+1.);
         // CALCULATION OF THE ELECTROMAGNETIC FIELDS
         *Em=*Em+MM*GTE[jlm(l,m)]-NM*GTM[jlm(l,m)]; 
         *Ez=*Ez+MZ*GTE[jlm(l,m)]-NZ*GTM[jlm(l,m)];
         *Ep=*Ep+MP*GTE[jlm(l,m)]-NP*GTM[jlm(l,m)];
         *Hm=*Hm+MM*GTM[jlm(l,m)]+NM*GTE[jlm(l,m)]; 
         *Hz=*Hz+MZ*GTM[jlm(l,m)]+NZ*GTE[jlm(l,m)];
         *Hp=*Hp+MP*GTM[jlm(l,m)]+NP*GTE[jlm(l,m)];
      }
   }
}
//------------------------------------------------------------------------------
// POSITION CALCULATIONS - LOOPS - COMPLETE CALCULATIONS
//------------------------------------------------------------------------------
void PVectorSphericalWaveFunctions(double *k,double *x, double *y, double *z,
                    int *lmax, int *nx, int *ny, int *nz, 
                    double complex *GTE, double complex *GTM,
                    double *rx, double *ry, double *rz,
                    double complex *Em, double complex *Ez, double complex *Ep,
                    double complex *Hm, double complex *Hz, double complex *Hp
                    ){ 
//------------------------------------------------------------------------------
   int i=0;
   for(int ix=0; ix<*nx; ix++){
      for(int iy=0; iy<*ny; iy++){
         for(int iz=0; iz<*nz; iz++){
           rx[i]=x[ix];
           ry[i]=y[iy];
           rz[i]=z[iz];
           VectorSphericalWaveFunctionsKernel(k,&x[ix],&y[iy],&z[iz],lmax,GTE,GTM,&Em[i],&Ez[i],&Ep[i],&Hm[i],&Hz[i],&Hp[i]);
           i++; 
         }    
      }    
   }    
//------------------------------------------------------------------------------
}
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
