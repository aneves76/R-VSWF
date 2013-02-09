#-------------------------------------------------------------------------------
Jnx.app<-function(n,x){
   return(((x/2)^n)/gamma(n+1))
}
Jnx.lim<-function(n){
   return(sqrt(n+1))
}
Jnx.def<-function(n,x){
   return(besselJ(x, n))
}   
#-------------------------------------------------------------------------------
jnx.app<-function(n,x){
   return(sqrt(pi/(2*x))*Jnx.app(n+.5,x))
}
jnx.lim<-function(n){
   return(sqrt(n/2))
}
jnx.def<-function(n,x){
   return(sqrt(pi/(2*x))*Jnx.def(n+.5,x))
}   
#-------------------------------------------------------------------------------
# LIMITE INVERTIDO
n.lim<-function(xm){
   return(2*xm^2)
}
#-------------------------------------------------------------------------------
#no<-3
#xo<-Jnx.lim(no)
#xj<-seq(0,xo,xo/200)
#ja<-jnx.app(no,xj)
#jd<-jnx.def(no,xj)
#plot(xj,ja,type='l')
#points(xj,jd,type='l',col='red')






