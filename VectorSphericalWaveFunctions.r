if(!is.loaded('VectorSphericalWaveFunctions'))
{
   cat("Loading shared library \"VectorSphericalWaveFunctions.so\".\n")
   dyn.load("VectorSphericalWaveFunctions.so")
#}else{
#   cat("Shared library \"VectorSphericalWaveFunctions.so\" already loaded.\n")
#   cat("Unloading shared library \"VectorSphericalWaveFunctions.so\".\n")
#   dyn.unload("VectorSphericalWaveFunctions.so")
#   cat("Loading shared library \"VectorSphericalWaveFunctions.so\".\n")
#   dyn.load("VectorSphericalWaveFunctions.so")
}
VectorSphericalWaveFunctions<-function(k,x,y,z,lmax,gte,gtm)
{
x[x==0]<-.Machine$double.xmin
y[y==0]<-.Machine$double.xmin
z[z==0]<-.Machine$double.xmin
.C("VectorSphericalWaveFunctions",
   k=as.double(k),
   x=as.double(x),
   y=as.double(y),
   z=as.double(z),
   lmax=as.integer(lmax),
   GTE=as.complex(gte),
   GTM=as.complex(gtm),
   Em=as.complex(0),
   Ez=as.complex(0),
   Ep=as.complex(0),
   Hm=as.complex(0),
   Hz=as.complex(0),
   Hp=as.complex(0)
   )
}
