#-------------------------------------------------------------------------------
# Testa recorrencia funcoes de Bessel Esférica
#-------------------------------------------------------------------------------
lmax<-300
x<-30
flm<-rep(1,lmax)
flp<-rep(1,lmax)
flm[1]<-0
flm[2]<--1/x
flp[1]<-1/x
flp[2]<-(1/x)*(3/x^2-1)
for(l in 2:lmax){
   flm[l+1]<-(2*l+1)*flm[l]/x-flm[l-1]
   flp[l+1]<-(2*l+1)*flp[l]/x-flp[l-1]
}
print(cbind(flm,flp))
