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
M<-7
N<-9
kx<-M*pi/a
ky<-N*pi/b
gama<-sqrt(kx^2+ky^2)
kz<-sqrt(k^2-gama^2)
dx<-a/(NP-1)
dy<-b/(NP-1)
x<-seq(0,a,by=dx)
y<-seq(0,b,by=dy)
z<-x
lmax<-as.integer(max(c(10,max(k*abs(x)),max(k*abs(y)))))
#lmax<-lmax*2
cat(SEP)
cat("LMAX=",lmax,"\n",sep="")
#-------------------------------------------------------------------------------
#xo<-sample(x,1)
#yo<-sample(y,1)
#zo<-sample(z,1)
xo<-x[111]+dx/2
yo<-y[171]+dy/2
zo<-z[182]+dx/2
#xo<-yo<-zo<-0
#-------------------------------------------------------------------------------
x<-x-xo
y<-y-yo
z<-z-zo
#-------------------------------------------------------------------------------
#xe<-sample(x,1)
#ye<-sample(y,1)
#ze<-sample(z,1)
xe<-x[38]
ye<-y[64]
ze<-z[15]
#-------------------------------------------------------------------------------
S<-1
P<-1
M<-4
#-------------------------------------------------------------------------------
Em<-((1+P)*psi.wg(gama,kz,xe+xo,ye+yo,ze+zo,(M-1),S)/2-P*((gama/k)^2)*psi.wg(gama,kz,xe+xo,ye+yo,ze+zo,(M-1)+S*(P-1),S)/2)
Ez<--1i*P*S*((gama/k)*(kz/k))*psi.wg(gama,kz,xe+xo,ye+yo,ze+zo,(M-1)+S*P,S)/sqrt(2)
Ep<-((1-P)*psi.wg(gama,kz,xe+xo,ye+yo,ze+zo,(M-1),S)/2+P*((gama/k)^2)*psi.wg(gama,kz,xe+xo,ye+yo,ze+zo,(M-1)+S*(P+1),S)/2)
#
Hm<--1i*P*(kz/k)*psi.wg(gama,kz,xe+xo,ye+yo,ze+zo,M-1,S)*(1+P)/2
Hz<--S*(gama/k)*psi.wg(gama,kz,xe+xo,ye+yo,ze+zo,(M-1)+S*P,S)/sqrt(2)
Hp<--1i*P*(kz/k)*psi.wg(gama,kz,xe+xo,ye+yo,ze+zo,M-1,S)*(1-P)/2
#-------------------------------------------------------------------------------
E<-c(Em,Ez,Ep)
H<-c(Hm,Hz,Hp)
cat(SEP)
print(cbind(E,H))
cat(SEP)
#-------------------------------------------------------------------------------
u<-BesselBeamP(gama,kz,xo,yo,zo,lmax,M,S,P)
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
cat("M =",M,"S =",S,"P =",P,"m = ",S*(M-1)-P,"\n")
cat(SEP)
print(data.frame(ReH=Re(H)/Re(H.pwe),ImH=Im(H)/Im(H.pwe),ReE=Re(E)/Re(E.pwe),ImE=Im(E)/Im(E.pwe)))
cat(SEP)
#-------------------------------------------------------------------------------
source("BesselBeamsP.r")
source("PVectorSphericalWaveFunctions.r")
#-------------------------------------------------------------------------------
# Guia de onda cilindrico
BBP<-BesselBeamsP(P,M,S,gama,kz,x+xo,y+yo,ze+zo)
#
tem.BBP<-array(BBP$Em,c(NP,NP,1))[,,1]
tez.BBP<-array(BBP$Ez,c(NP,NP,1))[,,1]
tep.BBP<-array(BBP$Ep,c(NP,NP,1))[,,1]
thm.BBP<-array(BBP$Hm,c(NP,NP,1))[,,1]
thz.BBP<-array(BBP$Hz,c(NP,NP,1))[,,1]
thp.BBP<-array(BBP$Hp,c(NP,NP,1))[,,1]
#
x11();image(x+xo,y+yo,Re(tem.BBP),main="BBP Em",col=cm.colors(1024));grid()
x11();image(x+xo,y+yo,Re(tez.BBP),main="BBP Ez",col=cm.colors(1024));grid()
x11();image(x+xo,y+yo,Re(tep.BBP),main="BBP Ep",col=cm.colors(1024));grid()
x11();image(x+xo,y+yo,Re(thm.BBP),main="BBP Hm",col=cm.colors(1024));grid()
x11();image(x+xo,y+yo,Re(thz.BBP),main="BBP Hz",col=cm.colors(1024));grid()
x11();image(x+xo,y+yo,Re(thp.BBP),main="BBP Hp",col=cm.colors(1024));grid()
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
