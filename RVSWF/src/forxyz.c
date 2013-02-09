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
