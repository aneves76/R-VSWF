v<-1:9
t<-c()
l<-c(1)
s<-1
for(i in 1:9){
   for(j in 1:9){
      for(k in 1:9){
         t<-c(t,i*100+j*10+k)
         l<-c(l,s)
         s<-s+1
      }
   }
}
T<-array(t,c(9,9,9))
get.ijk<-function(i,j,k){return(i*100+j*10+k)}
