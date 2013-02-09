#-------------------------------------------------------------------------------
# ESCREVER PARA CONVERTER LOCAL
rm(list=ls())
setwd("~/R-VSWF")
library(colorout)
SEP<-"#-------------------------------------------------------------------------------\n"
#-------------------------------------------------------------------------------
# FUNCOES
#-------------------------------------------------------------------------------
source("BSC.R")
source("HansenMultipoles.R")
#-------------------------------------------------------------------------------
ZJ <-matrix(c( 2.4048, 3.8317, 5.1356, 6.3802, 7.5883, 8.7715,
               5.5201, 7.0156, 8.4172, 9.7610,11.0647,12.3386,
               8.6537,10.1735,11.6198,13.0152,14.3725,15.7002,
              11.7915,13.3237,14.7960,16.2235,17.6160,18.9801,
              14.9309,16.4706,17.9598,19.4094,20.8269,22.2178),
     ncol=6,byrow=TRUE)
#-------------------------------------------------------------------------------
ZdJ<-matrix(c( 3.8317, 1.8412, 3.0542, 4.2012, 5.3175, 6.4156,
               7.0156, 5.3314, 6.7061, 8.0152, 9.2824,10.5199,
              10.1735, 8.5363, 9.9695,11.3459,12.6819,13.9872,
              13.3237,11.7060,13.1704,14.5858,15.9641,17.3128,
              16.4706,14.8636,16.3475,17.7887,19.1960,20.5755),
     ncol=6,byrow=TRUE)
#-------------------------------------------------------------------------------
# PARAMETROS DO GUIA DE ONDA
#-------------------------------------------------------------------------------
lambda=.5e-6
k=2*pi/lambda
n.lambda=10
NP=200
R<-5*lambda
M<-4
N<-2
#
x<-seq(-R,R,by=2*R/(NP-1))
y<-x
z<-x
radius<-function(x,y){return(sqrt(x^2+y^2))}
cosphi<-function(x,y){return(x/sqrt(x^2+y^2))}
sinphi<-function(x,y){return(y/sqrt(x^2+y^2))}
rho<-outer(x,y,FUN='radius')
cph<-outer(x,y,FUN='cosphi')
sph<-outer(x,y,FUN='sinphi')
phi.xy<-function(M,g,r,c,s){
   return(besselJ(M,g*r)*((c+1i*s)^M))
}
#-------------------------------------------------------------------------------
# MODOS: M <- M+1 já que J_0(x) ocupa a coluna 1
# TE - ZdJ
# TM - ZJ
#-------------------------------------------------------------------------------
g<-ZJ[N,M+1]/R
kz<-sqrt(k^2-g^2)
#U<-phi.xy(M,g,rho,cph,sph)
#image(x,y,Re(U))
##
lmax<-as.integer(max(c(10,max(k*R))))
lmax<-lmax*2
cat(SEP)
cat("LMAX=",lmax,"\n",sep="")
##-------------------------------------------------------------------------------
## POSICAO DA EXPANSAO (DE 1 A 200)
##-------------------------------------------------------------------------------
#nro<-sample(1:200,1)
#npo<-sample(1:200,1)
#nzo<-sample(1:200,1)
##xo<-r[nro]*cos(ph[npo])
##yo<-r[nro]*sin(ph[npo])
##zo<-z[nzo]
xo<-0
yo<-0
zo<-0
##-------------------------------------------------------------------------------
## POSICAO DO CAMPO TESTADO (DE 1 A 200)
##-------------------------------------------------------------------------------
xe<-sample(x,1)
ye<-sample(y,1)
ze<-sample(z,1)
if(sqrt(xe^2+ye^2+ze^2)>R){print("::::::::::::::::::::::::::::: r>R :::::::::::::::::::::::::::::")}
##-------------------------------------------------------------------------------
## Calcula os campos
##-------------------------------------------------------------------------------
## MODO TRANSVERSAL MAGNETICO s=1
# psi.wg(gama,kz,x,y,z,m,s=1)
E.TM.H.TE.m<-psi.wg(g,kz,xe+xo,ye+yo,ze+zo,M+1,s=1)*(-1i*kz/g)/sqrt(2)
E.TM.H.TE.z<-psi.wg(g,kz,xe+xo,ye+yo,ze+zo,M,s=1)
E.TM.H.TE.p<-psi.wg(g,kz,xe+xo,ye+yo,ze+zo,M-1,s=1)*( 1i*kz/g)/sqrt(2)
#
H.TM.E.TE.m<-psi.wg(g,kz,xe+xo,ye+yo,ze+zo,M+1,s=1)*(k/g)/sqrt(2)
H.TM.E.TE.z<-0
H.TM.E.TE.p<-psi.wg(g,kz,xe+xo,ye+yo,ze+zo,M-1,s=1)*(k/g)/sqrt(2)
#
E.TM.H.TE<-c(E.TM.H.TE.m,E.TM.H.TE.z,E.TM.H.TE.p)
H.TM.E.TE<-c(H.TM.E.TE.m,H.TM.E.TE.z,H.TM.E.TE.p)
FE<-sum(Conj(E.TM.H.TE)*E.TM.H.TE)
FH<-sum(Conj(H.TM.E.TE)*H.TM.E.TE)
#-------------------------------------------------------------------------------
cat(SEP)
print(matrix(c(k,g,kz,xo,yo,zo,xe,ye,ze,abs(xe-xo)/R,abs(ye-yo)/R,abs(ze-zo)/R),nrow=3)) 
cat(SEP)
#-------------------------------------------------------------------------------
cat(SEP)
print(cbind(E.TM.H.TE,H.TM.E.TE))
cat(SEP)
#-------------------------------------------------------------------------------
# Calcula os BSC
#-------------------------------------------------------------------------------
#MODE = 3, TE
#MODE = 4, TM
cat(SEP)
#WaveGuides(kx,ky,kz,xo,yo,zo,lmax,mode=0)
#WaveGuides(M ,N ,kz ,xo,yo,zo,lmax,mode=4)
u<-WaveGuides(M,N,kz,xo,yo,zo,lmax,mode=3)
# h<-rep(0,jlm(lmax,lmax))
# for(l in 1:lmax){h[jlm(l,M)]<-2*pi*(-1i)^M}
# u$GTE<-h*u$A
# u$GTM<-h*u$B
v<-HansenMultipoles(k,xe,ye,ze,lmax)
#
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
##-------------------------------------------------------------------------------
## GRAFICO
##-------------------------------------------------------------------------------
##cc<-function(x,y){
##   return(cos(kx*x)*cos(ky*y))
##}
##ccxy<-outer(x+xo,y+yo,FUN='cc')
##image(x,y,ccxy)
##grid()
##abline(v=xe,h=ye)
