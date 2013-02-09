if(!is.loaded('RectangularWaveGuide'))
{
   cat("Loading shared library \"RectangularWaveGuide.so\".\n")
   dyn.load("RectangularWaveGuide.so")
#}else{
#   cat("Shared library \"RectangularWaveGuide.so\" already loaded.\n")
#   cat("Unloading shared library \"RectangularWaveGuide.so\".\n")
#   dyn.unload("RectangularWaveGuide.so")
#   cat("Loading shared library \"RectangularWaveGuide.so\".\n")
#   dyn.load("RectangularWaveGuide.so")
}
RectangularWaveGuide<-function(TE=TRUE,kx,ky,kz,x,y,z){
if(TE)
{
   te<-1
}else{
   te<-0
}
nx<-length(x)
ny<-length(y)
nz<-length(z)
dummy<-rep(0,nx*ny*nz)
.C("RectangularWaveGuide",
   TE=as.integer(te),
   nx=as.integer(nx),ny=as.integer(ny),nz=as.integer(nz),
   kx=as.double(kx) ,ky=as.double(ky) ,kz=as.double(kz),
   x=as.double(x) ,y=as.double(y) ,z=as.double(z),
   rx=as.double(dummy) ,ry=as.double(dummy) ,rz=as.double(dummy),
   Hm=as.complex(dummy),Hz=as.complex(dummy),Hp=as.complex(dummy),
   Em=as.complex(dummy),Ez=as.complex(dummy),Ep=as.complex(dummy)
   )
}
#-------------------------------------------------------------------------------
# int *TE,
# int *nx, int *ny, int *nz,
# double *kx, double *ky, double *kz,
# double *x,  double *y,  double *z,
# double *rx, double *ry, double *rz,
# double complex *Hm, double complex *Hz, double complex *Hp,
# double complex *Em, double complex *Ez, double complex *Ep
#-------------------------------------------------------------------------------
