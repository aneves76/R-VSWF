GetFields<-function(A,nx,ny,nz){
   EM<-array(A$Em,c(nx,ny,nz))
   EZ<-array(A$Ez,c(nx,ny,nz))
   EP<-array(A$Ep,c(nx,ny,nz))
   HM<-array(A$Hm,c(nx,ny,nz))
   HZ<-array(A$Hz,c(nx,ny,nz))
   HP<-array(A$Hp,c(nx,ny,nz))
}
