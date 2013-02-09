#-------------------------------------------------------------------------------
# ESCREVER QUALQUER COISA
# INDICES
#-------------------------------------------------------------------------------
jlm<-function(l,m){
   j.lm<-l*(l+1)+m+1    
   j.tf<-abs(m)>l
   j.lm[j.tf]<-1
   return(j.lm)
}
#-------------------------------------------------------------------------------
# POLINOMIOS DE LEGENDRE E DERIVADA
#-------------------------------------------------------------------------------
LegendrePolynomialsQlm<-function(cth,lmax){
   sth<-sqrt(1-cth^2)
#-------------------------------------------------------------------------------
   LMAX=lmax*(lmax+2)+1
#-------------------------------------------------------------------------------
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
   dQlm[jlm(1, 1)]<--(cth/sth^2)*Qlm[jlm(1, 1)]                           # dQ11
   dQlm[jlm(1, 0)]<--(cth*Qlm[jlm(1, 0)]-sqrt(3)*Qlm[jlm(0,0)])/(sth^2)   # dQ10
   dQlm[jlm(1,-1)]<--(cth/sth^2)*Qlm[jlm(1,-1)]                           # dQ11
#-------------------------------------------------------------------------------
   if(lmax>1){
      for(l in 2:lmax){
         m.p<-0:(l-2)
         m.m<-(-l):(-1)
         mm<-(-l):l
         csc2<-1/(sth^2)
         csp<-(-1)^m.m
         Qlm[jlm(l,l  )]=-gammaQ(l)*sth*Qlm[jlm(l-1,l-1)]                                  #OK
         Qlm[jlm(l,l-1)]=deltaQ(l)*cth*Qlm[jlm(l-1,l-1)]                                   #OK
         Qlm[jlm(l,m.p)]=alfaQ(l,m.p)*cth*Qlm[jlm(l-1,m.p)]-betaQ(l,m.p)*Qlm[jlm(l-2,m.p)] #OK 
         Qlm[jlm(l,m.m)]=csp*Qlm[jlm(l,abs(m.m))]                                          #OK
         dQlm[jlm(l,mm)]<--l*cth*csc2*Qlm[jlm(l,mm)]+sqrt((2*l+1)/(2*l-1))*sqrt((l+mm)*(l-mm))*csc2*Qlm[jlm(l-1,mm)]
      }
   } 
   return(data.frame(Qlm,dQlm))
}
#-------------------------------------------------------------------------------
# ONDA PLANA 
#-------------------------------------------------------------------------------
PlaneWave<-function(kx,ky,kz,ux,uy,uz,lmax,xo=0,yo=0,zo=0){
   a<-SphericalHarmonics(kx,ky,kz,lmax) # Calculo direto de X_{lm}
   k<-sqrt(Conj(kx)*kx+Conj(ky)*ky+Conj(kz)*kz) # k real
   u<-sqrt(Conj(ux)*ux+Conj(uy)*uy+Conj(uz)*uz) # u imag
   LMAX=lmax*(lmax+2)+1
   # NORMALIZACAO k
   hkx<-kx/k
   hky<-ky/k
   hkp<-(hkx+1i*hky)/sqrt(2)
   hkm<-(hkx-1i*hky)/sqrt(2)
   hkz<-kz/k
   # NORMALIZACAO u
   hux<-ux/u
   huy<-uy/u
   huz<-uz/u
   hup<-(hux+1i*huy)/sqrt(2)
   hum<-(hux-1i*huy)/sqrt(2)
   # ALFA a = k x u
   ham<-1i*(hkm*huz-hkz*hum)
   haz<-1i*(hkp*hum-hkm*hup)
   hap<-1i*(hkz*hup-hkp*huz)
   # POLINOMIOS DE LEGENDRE
   U<-LegendrePolynomialsQlm(kz/k,lmax+1)
   qlm<-U$Qlm
   rm(U)
   # CONSTANTES
   cplm<-function(l,m){return(sqrt(l*(l+1)-m*(m+1)))}
   cmlm<-function(l,m){return(sqrt(l*(l+1)-m*(m-1)))}
   il<-function(l){return((1i^l)/sqrt(l*(l+1)))}
   # CALCULA
   GTE<-rep(0,LMAX)
   GTM<-rep(0,LMAX)
   Il<-rep(0,LMAX)
   for(l in 1:lmax){
      m<--l:l
      lmm<-jlm(l,m-1)
      lmz<-jlm(l,m)
      lmp<-jlm(l,m+1)
      Il[lmz]<-1i^l
      if((kx==0)&&(ky==0)){
         em<-1
         ez<-1
         ep<-1
      }else{
         gk<-sqrt(kx^2+ky^2)
         czt<-kx/gk
         szt<-ky/gk
         em<-(czt-1i*szt)^(m-1)
         ez<-(czt-1i*szt)^m
         ep<-(czt-1i*szt)^(m+1)
      }
      GTE[lmz]<-4*pi*il(l)*(cplm(l,m)*ep*qlm[lmp]*hup/sqrt(2)+m*ez*qlm[lmz]*huz+cmlm(l,m)*em*qlm[lmm]*hum/sqrt(2))
      GTM[lmz]<-4*pi*il(l)*(cplm(l,m)*ep*qlm[lmp]*hap/sqrt(2)+m*ez*qlm[lmz]*haz+cmlm(l,m)*em*qlm[lmm]*ham/sqrt(2))
   } 
   GTE0<-4*pi*Il*(Conj(a$Xm)*hum+Conj(a$Xz)*huz+Conj(a$Xp)*hup)
   GTM0<-4*pi*Il*(Conj(a$Xm)*ham+Conj(a$Xz)*haz+Conj(a$Xp)*hap)
   return(data.frame(GTE,GTM,GTE0,GTM0))
}
#-------------------------------------------------------------------------------
# ONDA PLANA NA DIRECAO z
#-------------------------------------------------------------------------------
PlaneWave.z<-function(lmax,norm=FALSE,s=1){
#-------------------------------------------------------------------------------
   if(lmax<1){lmax<-1}                         # Pelo menos 1 termo
   LMAX=lmax*(lmax+2)+1                        # Vetor para lmax
#-------------------------------------------------------------------------------
   gte<-rep(0,LMAX)
   gtm<-rep(0,LMAX)
   glm<-function(l,m,norm=FALSE){
      k<-(1i^l)*sqrt(4*pi*(2*l+1))
      k<-rep(k,2*l+1)
      k[m!=s]<-0
      if(norm){
         k<-k/sqrt(2)
      }
      return(k)
   }
   for(l in 0:lmax){
      m<--l:l
      gte[jlm(l,m)]<-glm(l,m,norm)
   }
   gtm<--1i*s*gte
   U<-data.frame(GTE=gte,GTM=gtm)
   return(U)
}
#-------------------------------------------------------------------------------
# GUIAS DE ONDA
#-------------------------------------------------------------------------------
WaveGuides<-function(kx,ky,kz,xo,yo,zo,lmax,mode=0,s=1){
   gama<-sqrt(kx^2+ky^2)
   k<-sqrt(gama^2+kz^2)
   LMAX=lmax*(lmax+2)+1
   #----------------------------------------
   u<-LegendrePolynomialsQlm(kz/k,lmax)
   Qlm<-u$Qlm
   dQlm<-u$dQlm
   rm(u)
   #----------------------------------------
   A<-rep(0,LMAX)
   B<-rep(0,LMAX)
   for(l in 0:lmax){
      m<--l:l
      if(l>0){
         llp1<-1/sqrt(l*(l+1))
      }else{
         llp1<-0
      }
      A[jlm(l,m)]<-2*(1i^l)*((k/gama)^2)*Qlm[jlm(l,m)]*m*llp1
      B[jlm(l,m)]<-2*(1i^(l-1))*dQlm[jlm(l,m)]*llp1
   }
   #----------------------------------------
   # MODOS PRINCIPAIS
   #----------------------------------------
   if(mode==0){
      cat("WAVE GUIDE TEST MODE\n")
      return(data.frame(A,B))
   }
   #----------------------------------------
   # GUIA DE ONDA RETANGULAR
   #----------------------------------------
   # TM WG
   if(mode==1){
      cat("RECTANGULAR WAVE GUIDE - TM MODE\n")
      Gf<-RectangularWG(kx,ky,kz,xo,yo,zo,lmax,s=-1)
      GTE<- A*Gf
      GTM<--B*Gf
      return(data.frame(GTE,GTM))
   }
   # TE WG
   if(mode==2){
      cat("RECTANGULAR WAVE GUIDE - TE MODE\n")
      Gf<-RectangularWG(kx,ky,kz,xo,yo,zo,lmax,s=1)
      GTE<-B*Gf
      GTM<-A*Gf
      return(data.frame(GTE,GTM))
   }
   if(mode==3){
      cat("CYLINDRICAL WAVE GUIDE - TE MODE\n")
      #WaveGuides(kx,ky,kz,xo,yo,zo,lmax,mode=0)
      #WaveGuides(M ,N  ,R,xo,yo,zo,lmax,mode=3)
      #CylindricalWG(m,n,r,x,y,z,lmax,s=signal(R),TM=FALSE){
      Gf<-CylindricalWG(kx,ky,kz,xo,yo,zo,lmax,s)
      GTE<-B*Gf
      GTM<-A*Gf
      return(data.frame(GTE,GTM))
   }
   if(mode==4){
      cat("CYLINDRICAL WAVE GUIDE - TM MODE\n")
      #WaveGuides(kx,ky,kz,xo,yo,zo,lmax,mode=0)
      #WaveGuides(M ,N ,R ,xo,yo,zo,lmax,mode=4)
      #CylindricalWG(m,n,r,x,y,z,lmax,s=signal(R),TM=TRUE){
      Gf<-CylindricalWG(kx,ky,kz,xo,yo,zo,lmax,s)
      GTE<- A*Gf
      GTM<--B*Gf
      return(data.frame(GTE,GTM))
   }
}
#-------------------------------------------------------------------------------
# RECTANGULAR WAVE GUIDE
#-------------------------------------------------------------------------------
RectangularWG<-function(kx,ky,kz,x,y,z,lmax,s=1){
   #s=+1 -> TE MODE
   #s=-1 -> TM MODE
   LMAX=lmax*(lmax+2)+1
   g<-rep(0,LMAX)
   gama<-sqrt(kx^2+ky^2)
   EXmY<-exp(1i*(kx*x-ky*y))
   EXpY<-exp(1i*(kx*x+ky*y))
   czt<-kx/gama
   szt<-ky/gama
   for(l in 0:lmax){
      m<--l:l
      eimz<-(czt+1i*szt)^m    
      f<-(-1)^m
      g[jlm(l,m)]<-pi*exp(1i*kz*z)*((EXmY+f*Conj(EXmY))*eimz+s*((EXpY+f*Conj(EXpY))*Conj(eimz)))/2
   }
   return(g)
}
#-------------------------------------------------------------------------------
# CYLINDRICAL WAVE GUIDE
#-------------------------------------------------------------------------------
#WaveGuides(kx,ky,kz,xo,yo,zo,lmax,mode=0)
CylindricalWG<-function(kx,ky,kz,x,y,z,lmax,s){
   m<-abs(s)
   s<-sign(s)
   gama<-sqrt(kx^2+ky^2)
   LMAX=lmax*(lmax+2)+1
   g<-rep(0,LMAX)
   #----------------------------------------
   for(l in 1:lmax){
      q<--l:l
      g[jlm(l,q)]<-2*pi*((-1i*s)^q)*psi.wg(gama,kz,x,y,z,m-s*q,s)
   }
   return(g)
}
#-------------------------------------------------------------------------------
# PSI(k,r)
#-------------------------------------------------------------------------------
psi.wg<-function(gama,kz,x,y,z,m,s=1){
   rho<-sqrt(x^2+y^2)
   if((x==0)&&(y==0)){
      cph<-1
      sph<-0
   }else{
      cph<-x/rho
      sph<-y/rho   
   }
   u<-besselJ(gama*rho,m)*((cph+1i*s*sph)^m)*exp(1i*kz*z)
   return(u)
}
#-------------------------------------------------------------------------------
# BESSEL BEAM I
#-------------------------------------------------------------------------------
BesselBeamTE<-function(gama,kz,xo,yo,zo,lmax,M,s){
   k<-sqrt(gama^2+kz^2)
   LMAX=lmax*(lmax+2)+1
   #----------------------------------------
   u<-LegendrePolynomialsQlm(kz/k,lmax+1)
   Qlm<-u$Qlm
   dQlm<-u$dQlm
   rm(u)
   #----------------------------------------
   GTE<-rep(0,LMAX)
   GTM<-rep(0,LMAX)
   for(l in 1:lmax){
      m<--l:l
      psio<-psi.wg(gama,kz,xo,yo,zo,M-s*m,s)
      if(l>0){
         A0<-4*pi*(1i^(l-s*m))*psio/sqrt(l*(l+1))
      }else{
         A0<-0
      }
     GTE[jlm(l,m)]<-m*Qlm[jlm(l,m)]*A0
     GTM[jlm(l,m)]<-1i*((gama/k)^2)*dQlm[jlm(l,m)]*A0
   }
   return(data.frame(GTE,GTM))
}
#-------------------------------------------------------------------------------
# BESSEL BEAM II
#-------------------------------------------------------------------------------
BesselBeamP<-function(gama,kz,xo,yo,zo,lmax,M,s,p){
   k<-sqrt(gama^2+kz^2)
   cth<-kz/k
   sth<-gama/k
   LMAX=lmax*(lmax+2)+1
   #----------------------------------------
   u<-LegendrePolynomialsQlm(kz/k,lmax+1)
   Qlm<-u$Qlm
   dQlm<-u$dQlm
   rm(u)
   #----------------------------------------
   clmp<-function(l,m,p){return(sqrt(l*(l+1)-m*(m+p)))}
   Klmqp<-function(l,m,q,p){return(((l+p*m)*(l+p*q))/((2*l-1)*(2*l+1)))}
   #----------------------------------------
   GTE<-rep(0,LMAX)
   GTM<-rep(0,LMAX)
   for(l in 1:lmax){
      m<--l:l
      psio<-psi.wg(gama,kz,xo,yo,zo,M-1-s*(m-p),s)
      if(l>0){
         A0<-4*pi*(1i^(l-s*(m-p)))*psio/sqrt(2*l*(l+1))
      }else{
         A0<-0
      }
     GTE[jlm(l,m)]<-A0*clmp(l,m,-p)*Qlm[jlm(l,m-p)]
     GTM[jlm(l,m)]<--1i*A0*(cth*sth*dQlm[jlm(l,m)]-P*m*Qlm[jlm(l,m)]/sth)  
   }
   return(data.frame(GTE,GTM))
}
#-------------------------------------------------------------------------------
# FIM
#-------------------------------------------------------------------------------
