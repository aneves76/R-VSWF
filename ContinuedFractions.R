#-------------------------------------------------------------------------------
# CALCULO DAS FUNCOES ESFERICAS DE BESSEL VIA CONTINUED FRACTIONS
# LENTZ METHOD
# fn=bo+(a1/b1+)(a2/b2+)(a3/b3+)...(an/bn)
#-------------------------------------------------------------------------------
#rm(list=ls())
#-------------------------------------------------------------------------------
eo<-.Machine$double.xmin
Tk<-function(x,k){return((2*k+1)/x)}
ACC<-10^-50
#-------------------------------------------------------------------------------
# L > x
# Opcao para x>L: Calcular jlm para um L2=int(1.5x)
#-------------------------------------------------------------------------------
NN<-50000
xo<-45000
print(NN/xo)
#-------------------------------------------------------------------------------
fn<-NN/xo
if(fn==0){fn<-eo}
Cn<-fn
Dn<-0
N<-1
DN<-10
fna<-c(fn)
while(abs(DN-1)>ACC){
   an<--1
   bn<-Tk(xo,N)
   Cn<-bn+an/Cn
   if(Cn==0){Cn<-eo}
   Dn<-bn+an*Dn
   if(Dn==0){Dn<-eo}
   Dn<-1/Dn
   DN<-Cn*Dn
   fn<-fn*DN
   N<-N+1
   fna<-c(fn,fna)
   if(N>2000){
      break
   }
}
print(N)
#-------------------------------------------------------------------------------
# SPHERICAL BESSEL FUNCTIONS CALCULATION
# DOWNWARD RECURRENCE 
# VALIDO PARA lmax>x
f<-fn
gn<-c(1)
dgn<-c(f)
gno<-1
dgno<-f  # dgno=gno*dgno
RN<-1
for(n in NN:1){
   gnom<-gno*(n+1)/xo+dgno
   dgno<-gnom*(n-1)/xo-gno
   gno<-gnom
   gn<-c(gno,gn)
   dgn<-c(dgno,dgn)
   if(abs(gno)>1e100){
      cat("RENORMALIZACAO ",RN,"N =",n,"\n")
      RN<-RN+1
      gn<-gn/gno
      dgn<-dgn/gno
      dgno<-dgno/gno
      gno<-1
   }
}
#-------------------------------------------------------------------------------
# NORMALIZACAO
gn<-(sin(xo)/xo)*gn/gn[1]
dgn<--dgn*gn[2]/dgn[1]
uo<-bessel_jl_steed_array(lmax=NN,x=xo)
if(!(TRUE%in%is.nan(gn))||(!TRUE%in%is.nan(uo))){
   plot(uo,pch=3,col='red',type='b')
   points(gn,type='b')
}else{
   cat("PROBLEMS\n")
   plot(fna,type='b')
}
#-------------------------------------------------------------------------------
# Comparacao da renormalizacao para a derivada
# a derivada e calculada por f_n'=(n/x)f_n-f_{n+1}
nox<-seq(0,NN)/xo 
dg<-(nox*gn)[1:NN]-gn[2:(NN+1)]
plot(dg,type='b')
points(dgn,pch=3,col='red')
