if(!is.loaded('PVectorSphericalWaveFunctions'))
{
   cat("Loading shared library \"PVectorSphericalWaveFunctions.so\".\n")
   dyn.load("PVectorSphericalWaveFunctions.so")
#}else{
#   cat("Shared library \"PVectorSphericalWaveFunctions.so\" already loaded.\n")
#   cat("Unloading shared library \"PVectorSphericalWaveFunctions.so\".\n")
#   dyn.unload("PVectorSphericalWaveFunctions.so")
#   cat("Loading shared library \"PVectorSphericalWaveFunctions.so\".\n")
#   dyn.load("PVectorSphericalWaveFunctions.so")
}
PVectorSphericalWaveFunctions<-function(k,x,y,z,lmax,gte,gtm)
{
x[x==0]<-.Machine$double.xmin
y[y==0]<-.Machine$double.xmin
z[z==0]<-.Machine$double.xmin
nx<-length(x)
ny<-length(y)
nz<-length(z)
dummy<-rep(0,nx*ny*nz)
.C("PVectorSphericalWaveFunctions",
   k=as.double(k),
   x=as.double(x),
   y=as.double(y),
   z=as.double(z),
   lmax=as.integer(lmax),
   nx=as.integer(nx),
   ny=as.integer(ny),
   nz=as.integer(nz),
   GTE=as.complex(gte),
   GTM=as.complex(gtm),
   rx=as.double(dummy),
   ry=as.double(dummy),
   rz=as.double(dummy),
   Em=as.complex(dummy),
   Ez=as.complex(dummy),
   Ep=as.complex(dummy),
   Hm=as.complex(dummy),
   Hz=as.complex(dummy),
   Hp=as.complex(dummy)
   )
}
