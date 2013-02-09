rm(list=ls())
source("PlaneWave.R")
SEP<-"================================================================================\n"
cat(SEP)
source("LegendrePolynomials.R")
cat(SEP)
rho<-sqrt(xo^2+yo^2)
r<-sqrt(rho^2+zo^2)
cth<-zo/r
sth<-rho/r
cph<-xo/rho
sph<-yo/rho
emph<-c(1)
for(l in 1:lmax){
   m<--l:l
   emph<-c(emph,(cph+1i*sph)^m)
}
cat(SEP)
print(c(cth,sth))
cat(SEP)
source("VectorSphericalWaveFunctions.r")
cat(SEP)
UVS<-VectorSphericalWaveFunctions(k,xo,yo,zo,lmax,u1$GTE,u1$GTM)
cat(SEP)
# u<-LegendrePolynomialsQlm(cth,lmax)
# u<-cbind(u,Eimph=emph,Ylm=u$Qlm*emph)
# u<-u[2:(lmax*(lmax+2)+1),]
# print(u)
E<-c(UVS$Em,UVS$Ez,UVS$Ep)
H<-c(UVS$Hm,UVS$Hz,UVS$Hp)
#print(cbind(v$N.m,v$N.z,v$N.p))
print(data.frame(E,H))
S<-SphericalHarmonics(xo,yo,zo,lmax)
cat(SEP)
