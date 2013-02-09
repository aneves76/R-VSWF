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
// CONSTANTES PARA Y E V
//------------------------------------------------------------------------------
double Km(int lo, int mo, int qo){
   double l=lo;
   double m=mo;
   double q=qo;
   return sqrt((l-m)*(l-q))/sqrt((2*l-1)*(2*l+1));
}
//------------------------------------------------------------------------------
double Ko(int lo, int mo, int qo){
   double l=lo;
   double m=mo;
   double q=qo;
   return sqrt((l-m)*(l+q))/sqrt((2*l-1)*(2*l+1));
}
//------------------------------------------------------------------------------
double Kp(int lo, int mo, int qo){
   double l=lo;
   double m=mo;
   double q=qo;
   return sqrt((l+m)*(l+q))/sqrt((2*l-1)*(2*l+1));
}
//------------------------------------------------------------------------------
//void CKonst(int lmax, double *cma, double *cpa, double *cza,
//      double *KmLm, double *KoLo, double *KpLM,
//      double *KmlM, double *Kolo, double *Kplm){
//   int i=0;
//   int l, m;
//   for(l=0;l<=lmax;l++){
//      for(m=-l;m<=l;m++){
//         cma[i]=cm(l,m);
//         cza[i]=m;
//         cpa[i]=cp(l,m);
//         KmLm[i]=Km(l+1,m,m-1);
//         KoLo[i]=Ko(l+1,m,m  );
//         KpLM[i]=Kp(l+1,m,m+1);
//         KmlM[i]=Km(l  ,m,m+1);
//         Kolo[i]=Ko(l  ,m,m  );
//         Kplm[i]=Kp(l  ,m,m-1);
//         i++;
//      }
//   }
//}
////------------------------------------------------------------------------------
//// CALCULATION OF CONSTANTS - ONE MORE LOOP - POSITION INDEPENDENT
////------------------------------------------------------------------------------
//// L DEPENDENT
//double CL[LMAX];
//double LL[LMAX];
//// CONSTANTS FOR X
//double CM[LMAX];
//double CZ[LMAX];
//double CP[LMAX];
//// CONSTANTS FOR Y AND V
//KmLm[LMAX];
//KoLo[LMAX];
//KpLM[LMAX];
//KmlM[LMAX];
//Kolo[LMAX];
//Kplm[LMAX];
//// SINGLE LOOP
//for(int l=1;l<=lmax;l++){
//   for(int m=-l; m<=l; m++){
//      CL[jlm(l,m)]=l;
//      LL[jlm(l,m)]=2*l+1;
//      CM[jlm(l,m)]=cm(l,m);
//      CZ[jlm(l,m)]=m;
//      CP[jlm(l,m)]=cp(l,m);
//      KmLm[jlm(l,m)]=Km(l+1,m,m-1);
//      KoLo[jlm(l,m)]=Ko(l+1,m,m  );
//      KpLM[jlm(l,m)]=Kp(l+1,m,m+1);
//      KmlM[jlm(l,m)]=Km(l  ,m,m+1);
//      Kolo[jlm(l,m)]=Ko(l  ,m,m  );
//      Kplm[jlm(l,m)]=Kp(l  ,m,m-1);
//   }
//}
//------------------------------------------------------------------------------
// POSITION DEPENDENT CALCULATIONS - MULTIPLES LOOPS - COMPLETE CALCULATIONS
//------------------------------------------------------------------------------
void VectorSphericalWaveFunctions(double *k,double *x, double *y, double *z,int *lmax, 
                    double complex *GTE, double complex *GTM,
                    double complex *Em, double complex *Ez, double complex *Ep,
                    double complex *Hm, double complex *Hz, double complex *Hp
                    ){
//   printf("%d\t%E\t%E\t%E\t%E\n",*lmax+1,*k,*x,*y,*z);
   int l,m;
   int LMAX=*lmax*(*lmax+2);
   int LMAXE=(*lmax+1)*(*lmax+3);
   double cph; 
   double sph;
   double rho=sqrt(*x*(*x)+*y*(*y));
   double r=sqrt(rho*rho+*z*(*z));
   double sth=rho/r;
   double cth=*z/r;
   if((*x==0)&&(*y==0)){
      cph=1;
      sph=0;
   }else{
      cph=*x/rho;
      sph=*y/rho;
   }
   // Spherical Bessel Funtions
   double JLM[*lmax+2];
  // double *JLM=&JLM0[0];
   gsl_sf_bessel_jl_steed_array(*lmax+1,*k*r,JLM);
   /* CALCULATIONS OK
   for(l=0;l<(*lmax+2);l++){
         printf("%d\t%f\t%E\n",l,*k*r,JLM[l]);
   } 
   */
   // Qlm - primeiros 4 termos
   double Qlm[LMAXE];
   Qlm[jlm(0, 0)]=1/sqrt(4*M_PI);
   Qlm[jlm(1, 1)]=-gammaQ(1)*sth*Qlm[jlm(0,0)]; // Q11
   Qlm[jlm(1, 0)]=sqrt(3.0)*cth*Qlm[jlm(0,0)];  // Q10
   Qlm[jlm(1,-1)]=-Qlm[jlm(1,1)];               // Q11*(-1)
   // Complex Exponencial for m=-1,0,1
   double complex Eim[2*(*lmax)+3];
   Eim[*lmax-1]=(cph-I*sph);
   Eim[*lmax  ]=1+I*0;
   Eim[*lmax+1]=(cph+I*sph);
   // Ylm - primeiros 4 termos
   double complex Ylm[LMAXE];
   Ylm[jlm(0, 0)]=Qlm[jlm(0, 0)];
   Ylm[jlm(1,-1)]=Qlm[jlm(1,-1)]*Eim[*lmax-1];
   Ylm[jlm(1, 0)]=Qlm[jlm(1, 0)];
   Ylm[jlm(1, 1)]=Qlm[jlm(1, 1)]*Eim[*lmax+1];
   /* OK jl, Qlm, Ylm
   for(l=0;l<2;l++){
      for(m=-l;m<=l;m++){
         printf("%d\t%d\t%d\t%f\t%f\t%f+%fi\n",l,m,jlm(l,m),JLM[jlm(l,m)],Qlm[jlm(l,m)],creal(Ylm[jlm(l,m)]),cimag(Ylm[jlm(l,m)]));
      }
   }
   printf("======================================================================\n");
   */
   // VECTOR SPHERICAL HARMONICS
   double complex XM; //[LMAX];
   double complex XZ; //[LMAX];
   double complex XP; //[LMAX];
   double complex YM; //[LMAX];
   double complex YZ; //[LMAX];
   double complex YP; //[LMAX];
   double complex VM; //[LMAX];
   double complex VZ; //[LMAX];
   double complex VP; //[LMAX];
   // HANSEN MULTIPOLES
   double complex MM; //[LMAX];
   double complex MZ; //[LMAX];
   double complex MP; //[LMAX];
   double complex NM; //[LMAX];
   double complex NZ; //[LMAX];
   double complex NP; //[LMAX];
   // OTHERS
   double kl;
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
      for(m=0; m<l; m++){
      // Demais valores de Qlm e Ylm
         Qlm[jlm(l+1, m)]=alfaQ(l+1,m)*cth*Qlm[jlm(l,m)]-betaQ(l+1,m)*Qlm[jlm(l-1,m)];
         Qlm[jlm(l+1,-m)]=pow(-1,m)*Qlm[jlm(l+1, m)];
         Ylm[jlm(l+1, m)]=Qlm[jlm(l+1, m)]*Eim[*lmax+m];
         Ylm[jlm(l+1,-m)]=Qlm[jlm(l+1,-m)]*Eim[*lmax-m];
      }
//------------------------------------------------------------------------------
//      for(m=-(l+1); m<=(l+1); m++){
//         Ylm[jlm(l+1,m)]=Qlm[jlm(l+1,m)]*Eim[*lmax+m];
//      }
//      for(m=-l; m<=l; m++){
//         printf("%d\t%d\t%d\t%f\t%f+%fi\t%f+%fi\n",l,m,jlm(l,m)+1,Qlm[jlm(l,m)],creal(Eim[*lmax+m]),cimag(Eim[*lmax+m]),creal(Ylm[jlm(l,m)]),cimag(Ylm[jlm(l,m)]));
//      }
//------------------------------------------------------------------------------
      for(int m=-l; m<=l; m++){
         XM=kl*cm(l,m)*Ylm[jlm(l,m-1)]/sqrt(2);
         XZ=kl*m*Ylm[jlm(l,m  )];
         XP=kl*cp(l,m)*Ylm[jlm(l,m+1)]/sqrt(2);
//         printf("--------X--------------------------------\n");
         printf("%d\t%d\t%d\t%f+%fi\n",l,m,jlm(l,m),creal(XM),cimag(XM));
//         printf("%d\t%d\t%d\t%f+%fi\n",l,m,jlm(l,m),creal(XZ),cimag(XZ));
//         printf("%d\t%d\t%d\t%f+%fi\n",l,m,jlm(l,m),creal(XP),cimag(XP));
         YM=(-Kp(l,m,m-1)*Ylm[jlm(l,m-1)]+Km(l+1,m,m-1)*Ylm[jlm(l+1,m-1)])/sqrt(2);
         YZ=  Ko(l,m,m  )*Ylm[jlm(l,m  )]+Ko(l+1,m,m  )*Ylm[jlm(l+1,m  )];
         YP=( Km(l,m,m+1)*Ylm[jlm(l,m+1)]-Kp(l+1,m,m+1)*Ylm[jlm(l+1,m+1)])/sqrt(2);
//        printf("--------Y--------------------------------\n");
//        printf("%d\t%d\t%d\t%f+%fi\n",l,m,jlm(l,m),creal(YM),cimag(YM));
//        printf("%d\t%d\t%d\t%f+%fi\n",l,m,jlm(l,m),creal(YZ),cimag(YZ));
//        printf("%d\t%d\t%d\t%f+%fi\n",l,m,jlm(l,m),creal(YP),cimag(YP));
         VM=kl*(-(l+1)*Kp(l,m,m-1)*Ylm[jlm(l,m-1)]-l*Km(l+1,m,m-1)*Ylm[jlm(l+1,m-1)])/sqrt(2);
         VZ=kl*( (l+1)*Ko(l,m,m  )*Ylm[jlm(l,m  )]-l*Ko(l+1,m,m  )*Ylm[jlm(l+1,m  )]);
         VP=kl*( (l+1)*Km(l,m,m+1)*Ylm[jlm(l,m+1)]+l*Kp(l+1,m,m+1)*Ylm[jlm(l+1,m+1)])/sqrt(2);
//        printf("--------V--------------------------------\n");
//        printf("%d\t%d\t%d\t%f+%fi\n",l,m,jlm(l,m),creal(VM),cimag(VM));
//        printf("%d\t%d\t%d\t%f+%fi\n",l,m,jlm(l,m),creal(VZ),cimag(VZ));
//        printf("%d\t%d\t%d\t%f+%fi\n",l,m,jlm(l,m),creal(VP),cimag(VP));
         // CALCULATION OF HANSEM MULTIPOLES
         MM=JLM[l]*XM;
         MZ=JLM[l]*XZ;
         MP=JLM[l]*XP;
//         printf("--------M--------------------------------\n");
//         printf("%d\t%d\t%d\t%f+%fi\n",l,m,jlm(l,m),creal(MM),cimag(MM));
//         printf("%d\t%d\t%d\t%f+%fi\n",l,m,jlm(l,m),creal(MZ),cimag(MZ));
//         printf("%d\t%d\t%d\t%f+%fi\n",l,m,jlm(l,m),creal(MP),cimag(MP));
         NM=((l+1)*JLM[l-1]-l*JLM[l+1])*VM/(2*l+1)+sqrt(l*(l+1))*JLM[l]*YM;
         NZ=((l+1)*JLM[l-1]-l*JLM[l+1])*VZ/(2*l+1)+sqrt(l*(l+1))*JLM[l]*YZ;
         NP=((l+1)*JLM[l-1]-l*JLM[l+1])*VP/(2*l+1)+sqrt(l*(l+1))*JLM[l]*YP;
//         printf("--------N--------------------------------\n");
//         printf("%d\t%d\t%d\t%f+%fi\n",l,m,jlm(l,m),creal(NM),cimag(NM));
//         printf("%d\t%d\t%d\t%f+%fi\n",l,m,jlm(l,m),creal(NZ),cimag(NZ));
//         printf("%d\t%d\t%d\t%f+%fi\n",l,m,jlm(l,m),creal(NP),cimag(NP));
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
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
