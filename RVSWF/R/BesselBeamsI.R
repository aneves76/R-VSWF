#-------------------------------------------------------------------------------
# ESCREVER PARA CONVERTER LOCAL
rm(list=ls())
graphics.off()
SEP<-"#-------------------------------------------------------------------------------\n"
#-------------------------------------------------------------------------------
source("BSC.R")
source("HansenMultipoles.R")
#-------------------------------------------------------------------------------
# PARAMETROS DO FEIXE DE BESSEL
#-------------------------------------------------------------------------------
lambda=.5e-6
k=2*pi/lambda
n.lambda=10
NP=200
a<-5*lambda
b<-8*lambda
M<-3
N<-5
kx<-M*pi/a
ky<-N*pi/b
gama<-sqrt(kx^2+ky^2)
kz<-sqrt(k^2-gama^2)
x<-seq(0,a,by=a/(NP-1))
y<-seq(0,b,by=b/(NP-1))
z<-x
lmax<-as.integer(max(c(10,max(k*abs(x)),max(k*abs(y)))))
#lmax<-lmax*2
cat(SEP)
cat("LMAX=",lmax,"\n",sep="")
#-------------------------------------------------------------------------------
xo<-sample(x,1)
yo<-sample(y,1)
zo<-sample(z,1)
#xo<-yo<-zo<-0
#-------------------------------------------------------------------------------
x<-x-xo
y<-y-yo
z<-z-zo
#-------------------------------------------------------------------------------
xe<-sample(x,1)
ye<-sample(y,1)
ze<-sample(z,1)
#-------------------------------------------------------------------------------
S<--1
M<-1
#-------------------------------------------------------------------------------
Em<- 1i*S*((kz/k)*(gama/k))*psi.wg(gama,kz,xe+xo,ye+yo,ze+zo,M-S,S)/sqrt(2)
Ez<-((gama/k)^2)*psi.wg(gama,kz,xe+xo,ye+yo,ze+zo,M,S)
Ep<--1i*S*((kz/k)*(gama/k))*psi.wg(gama,kz,xe+xo,ye+yo,ze+zo,M+S,S)/sqrt(2)
#
Hm<-S*(gama/k)*psi.wg(gama,kz,xe+xo,ye+yo,ze+zo,M-S,S)/sqrt(2)
Hz<-0
Hp<-S*(gama/k)*psi.wg(gama,kz,xe+xo,ye+yo,ze+zo,M+S,S)/sqrt(2)
#-------------------------------------------------------------------------------
E<-c(Em,Ez,Ep)
H<-c(Hm,Hz,Hp)
cat(SEP)
print(cbind(E,H))
cat(SEP)
#-------------------------------------------------------------------------------
u<-BesselBeamTE(gama,kz,xo,yo,zo,lmax,M,S)
v<-HansenMultipoles(k,xe,ye,ze,lmax)

Em.pwe<-sum(u$GTE*v$M.m-u$GTM*v$N.m)
Ez.pwe<-sum(u$GTE*v$M.z-u$GTM*v$N.z)
Ep.pwe<-sum(u$GTE*v$M.p-u$GTM*v$N.p)

Hm.pwe<-sum(u$GTM*v$M.m+u$GTE*v$N.m)
Hz.pwe<-sum(u$GTM*v$M.z+u$GTE*v$N.z)
Hp.pwe<-sum(u$GTM*v$M.p+u$GTE*v$N.p)
#-------------------------------------------------------------------------------
# CONFERINDO
#-------------------------------------------------------------------------------
E.pwe<-c(Em.pwe,Ez.pwe,Ep.pwe)
H.pwe<-c(Hm.pwe,Hz.pwe,Hp.pwe)
cat(SEP)
print(cbind(E.pwe,H.pwe))
cat(SEP)
#-------------------------------------------------------------------------------
source("BesselBeamsZ.r")
source("PVectorSphericalWaveFunctions.r")
#-------------------------------------------------------------------------------
# Guia de onda cilindrico
BBZ<-BesselBeamsZ(TE=TRUE,M,S,gama,kz,x+xo,y+yo,ze+zo)
#
tem.BBZ<-array(BBZ$Em,c(NP,NP,1))[,,1]
tez.BBZ<-array(BBZ$Ez,c(NP,NP,1))[,,1]
tep.BBZ<-array(BBZ$Ep,c(NP,NP,1))[,,1]
thm.BBZ<-array(BBZ$Hm,c(NP,NP,1))[,,1]
thz.BBZ<-array(BBZ$Hz,c(NP,NP,1))[,,1]
thp.BBZ<-array(BBZ$Hp,c(NP,NP,1))[,,1]
#
x11();image(x+xo,y+yo,Re(tem.BBZ),main="BBZ Em",col=cm.colors(1024));grid()
x11();image(x+xo,y+yo,Re(tez.BBZ),main="BBZ Ez",col=cm.colors(1024));grid()
x11();image(x+xo,y+yo,Re(tep.BBZ),main="BBZ Ep",col=cm.colors(1024));grid()
x11();image(x+xo,y+yo,Re(thm.BBZ),main="BBZ Hm",col=cm.colors(1024));grid()
x11();image(x+xo,y+yo,Re(thz.BBZ),main="BBZ Hz",col=cm.colors(1024));grid()
x11();image(x+xo,y+yo,Re(thp.BBZ),main="BBZ Hp",col=cm.colors(1024));grid()
#-------------------------------------------------------------------------------
# Partial Wave Expansion
UVS<-PVectorSphericalWaveFunctions(k,x,y,ze,lmax,u$GTE,u$GTM)
#
tem.uvs<-array(UVS$Em,c(NP,NP,1))[,,1]
tez.uvs<-array(UVS$Ez,c(NP,NP,1))[,,1]
tep.uvs<-array(UVS$Ep,c(NP,NP,1))[,,1]
thm.uvs<-array(UVS$Hm,c(NP,NP,1))[,,1]
thz.uvs<-array(UVS$Hz,c(NP,NP,1))[,,1]
thp.uvs<-array(UVS$Hp,c(NP,NP,1))[,,1]
#
x11();image(x+xo,y+yo,Re(tem.uvs),main=paste("PWE Em",lmax),col=cm.colors(1024));grid()
x11();image(x+xo,y+yo,Re(tez.uvs),main=paste("PWE Ez",lmax),col=cm.colors(1024));grid()
x11();image(x+xo,y+yo,Re(tep.uvs),main=paste("PWE Ep",lmax),col=cm.colors(1024));grid()
x11();image(x+xo,y+yo,Re(thm.uvs),main=paste("PWE Hm",lmax),col=cm.colors(1024));grid()
x11();image(x+xo,y+yo,Re(thz.uvs),main=paste("PWE Hz",lmax),col=cm.colors(1024));grid()
x11();image(x+xo,y+yo,Re(thp.uvs),main=paste("PWE Hp",lmax),col=cm.colors(1024));grid()
#-------------------------------------------------------------------------------



