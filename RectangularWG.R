#-------------------------------------------------------------------------------
# ESCREVER PARA CONVERTER LOCAL
rm(list=ls())
graphics.off()
SEP<-"#-------------------------------------------------------------------------------\n"
#-------------------------------------------------------------------------------
# ONDA PLANA
#-------------------------------------------------------------------------------
# Campo Eletrico \vec{k}=k\hat{k}_z
# E=Eo exp(ik_zz)
#-------------------------------------------------------------------------------
# FUNCOES
#-------------------------------------------------------------------------------
source("BSC.R")
source("HansenMultipoles.R")
#-------------------------------------------------------------------------------
# PARAMETROS DO GUIA DE ONDA
#-------------------------------------------------------------------------------
lambda=.5e-6
k=2*pi/lambda
n.lambda=10
NP=200
a<-5*lambda
b<-8*lambda
M<-7
N<-6
kx<-M*pi/a
ky<-N*pi/b
gamma<-sqrt(kx^2+ky^2)
kz<-sqrt(k^2-gamma^2)
dx<-a/(NP-1)
dy<-b/(NP-1)
x<-seq(0,a,dx)
y<-seq(0,b,dy)
z<-x
lmax<-as.integer(max(c(10,max(k*abs(x)),max(k*abs(y)))))
print(lmax)
cat(SEP)
cat("LMAX=",lmax,"\n",sep="")
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# POSICAO DA EXPANSAO (DE 1 A 200)
#-------------------------------------------------------------------------------
nxo<-sample(1:200,1)
nyo<-sample(1:200,1)
nzo<-sample(1:200,1)
# insere desvio em (xo,yo,zo) para evitar NaN
xo<-x[nxo]+dx/2
yo<-y[nyo]+dy/2
zo<-z[nzo]+dx/2
#
xo<-x[27]+dx/2
yo<-y[82]+dy/2
zo<-z[29]+dx/2
#xo<-0
#yo<-0
#zo<-0
#-------------------------------------------------------------------------------
# MUDA O REFERENCIAL
#-------------------------------------------------------------------------------
x<-x-xo
y<-y-yo
z<-z-zo
#-------------------------------------------------------------------------------
# POSICAO DO CAMPO TESTADO (DE 1 A 200)
#-------------------------------------------------------------------------------
nx<-sample(1:200,1)
ny<-sample(1:200,1)
nz<-sample(1:200,1)
xe<-x[nx]
ye<-y[ny]
ze<-z[nz]
xe<-x[73]
ye<-y[80]
ze<-z[173]
#xe<-5e-7
#ye<-0
#ze<-xe
#-------------------------------------------------------------------------------
# Calcula os campos
#-------------------------------------------------------------------------------
# MODO TRANSVERSAL MAGNETICO
E.TM.x<-1i*((kz*kx)/gamma^2)*cos(kx*(xe+xo))*sin(ky*(ye+yo))*exp(1i*kz*(ze+zo))
E.TM.y<-1i*((kz*ky)/gamma^2)*sin(kx*(xe+xo))*cos(ky*(ye+yo))*exp(1i*kz*(ze+zo))
E.TM.z<-sin(kx*(xe+xo))*sin(ky*(ye+yo))*exp(1i*kz*(ze+zo))
#
H.TM.x<--1i*((k*ky)/gamma^2)*sin(kx*(xe+xo))*cos(ky*(ye+yo))*exp(1i*kz*(ze+zo))
H.TM.y<- 1i*((k*kx)/gamma^2)*cos(kx*(xe+xo))*sin(ky*(ye+yo))*exp(1i*kz*(ze+zo))
H.TM.z<- 0
#
E.TM.m<-(E.TM.x-1i*E.TM.y)/sqrt(2)
E.TM.p<-(E.TM.x+1i*E.TM.y)/sqrt(2)
E.TM<-c(E.TM.m,E.TM.z,E.TM.p)
#
H.TM.m<-(H.TM.x-1i*H.TM.y)/sqrt(2)
H.TM.p<-(H.TM.x+1i*H.TM.y)/sqrt(2)
H.TM<-c(H.TM.m,H.TM.z,H.TM.p)
#-------------------------------------------------------------------------------
# MODO TRANSVERSAL ELETRICO
#-------------------------------------------------------------------------------
H.TE.x<--1i*((kz*kx)/gamma^2)*sin(kx*(xe+xo))*cos(ky*(ye+yo))*exp(1i*kz*(ze+zo))
H.TE.y<--1i*((kz*ky)/gamma^2)*cos(kx*(xe+xo))*sin(ky*(ye+yo))*exp(1i*kz*(ze+zo))
H.TE.z<- cos(kx*(xe+xo))*cos(ky*(ye+yo))*exp(1i*kz*(ze+zo))
#
E.TE.x<--1i*((k*ky)/gamma^2)*cos(kx*(xe+xo))*sin(ky*(ye+yo))*exp(1i*kz*(ze+zo))
E.TE.y<- 1i*((k*kx)/gamma^2)*sin(kx*(xe+xo))*cos(ky*(ye+yo))*exp(1i*kz*(ze+zo))
E.TE.z<- 0
#
#
E.TE.m<-(E.TE.x-1i*E.TE.y)/sqrt(2)
E.TE.p<-(E.TE.x+1i*E.TE.y)/sqrt(2)
E.TE<-c(E.TE.m,E.TE.z,E.TE.p)
#
H.TE.m<-(H.TE.x-1i*H.TE.y)/sqrt(2)
H.TE.p<-(H.TE.x+1i*H.TE.y)/sqrt(2)
H.TE<-c(H.TE.m,H.TE.z,H.TE.p)
#-------------------------------------------------------------------------------
cat(SEP)
cat("TM MODE\n")
cat(SEP)
print(cbind(E.TM,H.TM))
cat(SEP)
cat("TE MODE\n")
cat(SEP)
print(cbind(E.TE,H.TE))
#-------------------------------------------------------------------------------
# Calcula os BSC
#-------------------------------------------------------------------------------
#MODE = 1, TM
#MODE = 2, TE
cat(SEP)
u<-WaveGuides(kx,ky,kz,xo,yo,zo,lmax,mode=2)
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
source("RectangularWaveGuide.r")
source("PVectorSphericalWaveFunctions.r")
#-------------------------------------------------------------------------------
# Guia de onda cilindrico
RWG<-RectangularWaveGuide(TE=TRUE,kx,ky,kz,x+xo,y+yo,ze+zo)
#
tem.rwg<-array(RWG$Em,c(NP,NP,1))[,,1]
tez.rwg<-array(RWG$Ez,c(NP,NP,1))[,,1]
tep.rwg<-array(RWG$Ep,c(NP,NP,1))[,,1]
thm.rwg<-array(RWG$Hm,c(NP,NP,1))[,,1]
thz.rwg<-array(RWG$Hz,c(NP,NP,1))[,,1]
thp.rwg<-array(RWG$Hp,c(NP,NP,1))[,,1]
#
x11();image(x+xo,y+yo,Re(tem.rwg),main="RWG Em",col=cm.colors(1024));grid()
x11();image(x+xo,y+yo,Re(tez.rwg),main="RWG Ez",col=cm.colors(1024));grid()
x11();image(x+xo,y+yo,Re(tep.rwg),main="RWG Ep",col=cm.colors(1024));grid()
x11();image(x+xo,y+yo,Re(thm.rwg),main="RWG Hm",col=cm.colors(1024));grid()
x11();image(x+xo,y+yo,Re(thz.rwg),main="RWG Hz",col=cm.colors(1024));grid()
x11();image(x+xo,y+yo,Re(thp.rwg),main="RWG Hp",col=cm.colors(1024));grid()
#-------------------------------------------------------------------------------
# Partial Wave Expansion
UVS<-PVectorSphericalWaveFunctions(k,x,y,zo,lmax,u$GTE,u$GTM)
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
