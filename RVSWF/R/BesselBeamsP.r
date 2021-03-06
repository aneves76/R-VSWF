#-------------------------------------------------------------------------------
if(!is.loaded('BesselBeamsP'))
{
   cat("Loading shared library \"BesselBeamsP.so\".\n")
   dyn.load("BesselBeamsP.so")
#}else{
#   cat("Shared library \"BesselBeamsP.so\" already loaded.\n")
#   cat("Unloading shared library \"BesselBeamsP.so\".\n")
#   dyn.unload("BesselBeamsP.so")
#   cat("Loading shared library \"BesselBeamsP.so\".\n")
#   dyn.load("BesselBeamsP.so")
}
#-------------------------------------------------------------------------------
BesselBeamsP<-function(P,M,S,g,kz,x,y,z){
nx<-length(x)
ny<-length(y)
nz<-length(z)
dummy<-rep(0,nx*ny*nz)
.C("BesselBeamsP",
   P=as.integer(P),M=as.integer(M),S=as.integer(S),
   nx=as.integer(nx),ny=as.integer(ny),nz=as.integer(nz),
   gama=as.double(g),kz=as.double(kz),
   x=as.double(x),y=as.double(y),z=as.double(z),
   rx=as.double(dummy) ,ry=as.double(dummy) ,rz=as.double(dummy),
   Hm=as.complex(dummy),Hz=as.complex(dummy),Hp=as.complex(dummy),
   Em=as.complex(dummy),Ez=as.complex(dummy),Ep=as.complex(dummy)
   )
}
#-------------------------------------------------------------------------------
