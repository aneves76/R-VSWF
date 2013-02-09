#-------------------------------------------------------------------------------
# CALCULA OS HARMONICOS ESFERICOS ESCALARES (SSH)
#-------------------------------------------------------------------------------
LegendrePolynomialsQlm<-function(cth,lmax){
   sth<-sqrt(1-cth^2)
#-------------------------------------------------------------------------------
   LMAX=lmax*(lmax+2)+1
#-------------------------------------------------------------------------------
   jlm<-function(l,m){
      j.lm<-l*(l+1)+m+1
      j.tf<-abs(m)>l
      j.lm[j.tf]<-1
      return(j.lm)
   }
   alfaQ<-function(l,m){
      return(sqrt(((2*l-1)*(2*l+1))/((l-m)*(l+m))))
   }
   betaQ<-function(l,m){
      return(sqrt((2*l+1)/(2*l-3))*sqrt(((l+m-1)*(l-m-1))/((l-m)*(l+m))))
   }
   gammaQ<-function(l){
      return(sqrt((2*l+1)/(2*l)))
   }
   deltaQ<-function(l){
      return(sqrt(2*l+1))
   }
#-------------------------------------------------------------------------------
   if(lmax<2){
      Qlm<-rep(0,4)
      dQlm<-rep(0,4)
   }else{
      Qlm<-rep(0,LMAX)
      dQlm<-rep(0,LMAX)
   }
#-------------------------------------------------------------------------------
   Qlm[jlm(0,0)]<-1/sqrt(4*pi)                 #Q00
   Qlm[jlm(1,1)]<--gammaQ(1)*sth*Qlm[jlm(0,0)] #Q11
   Qlm[jlm(1,0)]<-sqrt(3)*cth*Qlm[1]           #Q10
   Qlm[jlm(1,-1)]<--Qlm[jlm(1,1)]              #Q11*(-1)
#-------------------------------------------------------------------------------
   dQlm[jlm(0, 0)]<-0                                                     # dQ00
   dQlm[jlm(1, 1)]<-(cth/sth)*Qlm[jlm(1, 1)]                              # dQ11
   dQlm[jlm(1, 0)]<-(cth/sth)*Qlm[jlm(1, 0)]-(sqrt(3)/sth)*Qlm[jlm(0,0)]  # dQ10
   dQlm[jlm(1,-1)]<-(cth/sth)*Qlm[jlm(1,-1)]                              # dQ11
#-------------------------------------------------------------------------------
   if(lmax>1){
      for(l in 2:lmax){
         m.p<-0:(l-2)
         m.m<-(-l):(-1)
         mm<--l:l
         csp<-(-1)^m.m
         Qlm[jlm(l,l)]=-gammaQ(l)*sth*Qlm[jlm(l-1,l-1)]                                    #OK
         Qlm[jlm(l,l-1)]=deltaQ(l)*cth*Qlm[jlm(l-1,l-1)]                                   #OK
         Qlm[jlm(l,m.p)]=alfaQ(l,m.p)*cth*Qlm[jlm(l-1,m.p)]-betaQ(l,m.p)*Qlm[jlm(l-2,m.p)] #OK 
         Qlm[jlm(l,m.m)]=csp*Qlm[jlm(l,abs(m.m))]                                          #OK
         dQlm[jlm(l,-l):jlm(l,l)]<-l*cth*Qlm[jlm(l,-l):jlm(l,l)]-sqrt((2*l+1)/(2*l-1))*sqrt((l+mm)*(l-mm))*Qlm[jlm(l-1,mm)]
      }
   } 
   return(data.frame(Qlm,dQlm))
}
