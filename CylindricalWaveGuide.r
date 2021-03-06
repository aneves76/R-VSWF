#-------------------------------------------------------------------------------
if(!is.loaded('CylindricalWaveGuide'))
{
   cat("Loading shared library \"CylindricalWaveGuide.so\".\n")
   dyn.load("CylindricalWaveGuide.so")
#}else{
#   cat("Shared library \"CylindricalWaveGuide.so\" already loaded.\n")
#   cat("Unloading shared library \"CylindricalWaveGuide.so\".\n")
#   dyn.unload("CylindricalWaveGuide.so")
#   cat("Loading shared library \"CylindricalWaveGuide.so\".\n")
#   dyn.load("CylindricalWaveGuide.so")
}
#-------------------------------------------------------------------------------
CylindricalWaveGuide<-function(TE=TRUE,M,S,g,kz,x,y,z){
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
.C("CylindricalWaveGuide",
   TE=as.integer(te),M=as.integer(M),S=as.integer(S),
   nx=as.integer(nx),ny=as.integer(ny),nz=as.integer(nz),
   gama=as.double(g),kz=as.double(kz),
   x=as.double(x),y=as.double(y),z=as.double(z),
   rx=as.double(dummy) ,ry=as.double(dummy) ,rz=as.double(dummy),
   Hm=as.complex(dummy),Hz=as.complex(dummy),Hp=as.complex(dummy),
   Em=as.complex(dummy),Ez=as.complex(dummy),Ep=as.complex(dummy)
   )
}
#-------------------------------------------------------------------------------
