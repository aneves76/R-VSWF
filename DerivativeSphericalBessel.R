#-------------------------------------------------------------------------------
# Comparacao da derivada de um vetor de n funcoes esfericas de Bessel
#-------------------------------------------------------------------------------
dfn<-function(x,fn){
   nmax<-length(fn)
   nma<-nmax-1
   n<-0:(length(fn)-1)
   np<-n+1
   np2<-2*n+1
   dfn1<-(n/np2)[2:nma]*fn[1:(nma-1)]-(np/np2)[2:nma]*fn[3:nmax]
   dfn1<-c(0,dfn1,0)
   dfn2<-fn[1:nma]-np[2:nmax]*fn[2:nmax]/x
   dfn2<-c(0,dfn2)
   dfn3<-n[1:nma]*fn[1:nma]/x-fn[2:nmax]
   dfn3<-c(dfn3,0)
   return(data.frame(dfn1,dfn2,dfn3))
}
