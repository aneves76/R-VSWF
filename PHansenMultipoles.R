#-------------------------------------------------------------------------------
# CALCULA A EXPANSAO PARA UM CONJUNTO DE PONTOS
# INPUT: k, x, y, z, lmax, GTE, GTM
#-------------------------------------------------------------------------------
SphericalHarmonics<-function(x,y,z,lmax){
#-------------------------------------------------------------------------------
# UNIDADES BASICAS
#-------------------------------------------------------------------------------
   rho<-sqrt(x^2+y^2)  # rho
   r<-sqrt(rho^2+z^2)  # r
   cth<-z/r            # cos(theta)
   sth<-rho/r          # sin(theta)
   cph<-x/rho          # cos(phi)
   sph<-y/rho          # sin(phi)
   if((x==0)&&(y==0)){
      cph<-1
      sph<-0
   }
#-------------------------------------------------------------------------------
# CONSTANTES PARA Qlm
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
# POSICIONAMENTO GERAL DOS ELEMENTOS EM RELACAO A (l,m) DOS MULTIPOLOS 
#-------------------------------------------------------------------------------
   jlm<-function(l,m){
      j.lm<-l*(l+1)+m+1    
      j.tf<-abs(m)>l
      j.lm[j.tf]<-1
      return(j.lm)
   }
#-------------------------------------------------------------------------------
# Constantes para X
#-------------------------------------------------------------------------------
   cp<-function(l,m){
      return(sqrt(l*(l+1)-m*(m+1)))
   }
   cm<-function(l,m){
      return(sqrt(l*(l+1)-m*(m-1)))
   }
   cz<-function(l,m){
      return(m)
   }
   cl<-function(l,m){
      return(l)
   }
   ll<-function(l,m){
      return(1/sqrt(l*(l+1)))
   }
#-------------------------------------------------------------------------------
# Constantes para Y e V
#-------------------------------------------------------------------------------
   Km<-function(l,m,q){
      return(sqrt((l-m)*(l-q))/sqrt((2*l-1)*(2*l+1)))
   }
   Ko<-function(l,m,q){
      return(sqrt((l-m)*(l+q))/sqrt((2*l-1)*(2*l+1)))
   }
   Kp<-function(l,m,q){
      return(sqrt((l+m)*(l+q))/sqrt((2*l-1)*(2*l+1)))
   }
#-------------------------------------------------------------------------------
# Constantes para Y.lm e Y.LM
#-------------------------------------------------------------------------------
   KYmlm<-function(l,m){
      return(sqrt((l+m-1)*(l+m))/sqrt(2*l*(2*l-1)))
   }
   KYolm<-function(l,m){
      return(sqrt((l-m)*(l+m))/sqrt(l*(2*l-1)))
   }
   KYplm<-function(l,m){
      return(sqrt((l-m-1)*(l-m))/sqrt(2*l*(2*l-1)))
   }
   #
   KYmLm<-function(l,m){
      return(sqrt((l-m+1)*(l-m+2))/sqrt(2*(l+1)*(2*l+3)))
   }
   KYoLm<-function(l,m){
      return(sqrt((l-m+1)*(l+m+1))/sqrt((l+1)*(2*l+3)))
   }
   KYpLm<-function(l,m){
      return(sqrt((l+m+2)*(l+m+1))/sqrt(2*(l+1)*(2*l+3)))
   }
#-------------------------------------------------------------------------------
# PRIMEIROS HARMONICOS ESFERICOS ESCALARES 
#-------------------------------------------------------------------------------
   if(lmax<1){lmax<-1}                         # Pelo menos 1 termo
   LMAX=lmax*(lmax+2)+1                        # Vetor para lmax
   LMAXE=(lmax+1)*(lmax+3)+1                   # Vetor para (lmax+1) - Extendido
#-------------------------------------------------------------------------------
   Qlm<-rep(0,LMAXE)                           # Primeiros 4 termos Qlm
   Qlm[jlm(0,0)]<-1/sqrt(4*pi)                 # Q00
   Qlm[jlm(1,1)]<--gammaQ(1)*sth*Qlm[jlm(0,0)] # Q11
   Qlm[jlm(1,0)]<-sqrt(3)*cth*Qlm[1]           # Q10
   Qlm[jlm(1,-1)]<--Qlm[jlm(1,1)]              # Q11*(-1)
#-------------------------------------------------------------------------------
   Ylm<-rep(1,LMAXE)
   eimphi.m<-rep(1,lmax+1)                 #m=0
   eimphi.p<-rep(1,lmax+1)
   eimphi.m[2]<-eimphi.m[1]*(cph-1i*sph)   #m=1
   eimphi.p[2]<-eimphi.p[1]*(cph+1i*sph)
#-------------------------------------------------------------------------------
   Ylm[1:4]<-Qlm[1:4]*c(1,(cph-1i*sph),1,(cph+1i*sph))
#-------------------------------------------------------------------------------
# VETORES 
#-------------------------------------------------------------------------------
   dummy<-rep(0,LMAX)
#------------------
   Y.l.m<-dummy
   Y.l.o<-dummy
   Y.l.M<-dummy
#------------------
   Y.o.m<-dummy
   Y.o.o<-dummy
   Y.o.M<-dummy
#------------------
   Y.L.m<-dummy
   Y.L.o<-dummy
   Y.L.M<-dummy
#------------------
# VECTOR SPH. HARM.
#------------------
   Vm<-dummy
   Vz<-dummy
   Vp<-dummy
#------------------
   Xm<-dummy
   Xz<-dummy
   Xp<-dummy
#------------------
   Ym<-dummy
   Yz<-dummy
   Yp<-dummy
#------------------
# CONSTANTES
#------------------
   c.l<-dummy
   l.l<-dummy
#------------------
   c.p<-dummy
   c.m<-dummy
   c.z<-dummy
#------------------
   KmLm<-dummy
   KoLo<-dummy
   KpLM<-dummy
#------------------
   KmlM<-dummy
   Kolo<-dummy
   Kplm<-dummy
#------------------
   Km.lm<-dummy
   Ko.lm<-dummy
   Kp.lm<-dummy
#------------------
   Km.Lm<-dummy
   Ko.Lm<-dummy
   Kp.Lm<-dummy
#-------------------------------------------------------------------------------
# NUCLEO ---  CALCULO DOS VETORES
#-------------------------------------------------------------------------------
   for(l in 1:lmax){
      mm<-(-l):l 
#-------------------------------------------------------------------------------
      # POLINOMIOS ASSOCIADOS DE LEGENDRE NORMALIZADOS E SSH
      l<-l+1                              #l inicia em 2 (deg 5)
      m.p<-0:(l-2)
      m.m<-(-l):(-1)
      csp<-(-1)^m.m
      Qlm[jlm(l,l  )]=-gammaQ(l)*sth*Qlm[jlm(l-1,l-1)]                                  #OK
      Qlm[jlm(l,l-1)]=deltaQ(l)*cth*Qlm[jlm(l-1,l-1)]                                   #OK
      Qlm[jlm(l,m.p)]=alfaQ(l,m.p)*cth*Qlm[jlm(l-1,m.p)]-betaQ(l,m.p)*Qlm[jlm(l-2,m.p)] #OK
      Qlm[jlm(l,m.m)]=csp*Qlm[jlm(l,abs(m.m))]                                          #OK
      # EXP(I*L*PHI)
      eimphi.p[l+1]<-eimphi.p[l]*(cph+1i*sph)
      eimphi.m[l+1]<-eimphi.m[l]*(cph-1i*sph)
      eimphi<-c(rev(eimphi.m[2:(l+1)]),eimphi.p[1:(l+1)])
      # Ylm
      Ylm[jlm(l,-l):jlm(l,l)]<-Qlm[jlm(l,-l):jlm(l,l)]*eimphi
      l<-l-1
#-------------------------------------------------------------------------------
# SSH DESLOCADOS
#-------------------------------------------------------------------------------
      Y.l.m[jlm(l,mm)]<-Ylm[jlm(l-1,mm-1)]
      Y.l.o[jlm(l,mm)]<-Ylm[jlm(l-1,mm  )]
      Y.l.M[jlm(l,mm)]<-Ylm[jlm(l-1,mm+1)]
#------------------------------------------------
      Y.o.m[jlm(l,mm)]<-Ylm[jlm(l  ,mm-1)]
      Y.o.o[jlm(l,mm)]<-Ylm[jlm(l  ,mm  )]
      Y.o.M[jlm(l,mm)]<-Ylm[jlm(l  ,mm+1)]
#------------------------------------------------
      Y.L.m[jlm(l,mm)]<-Ylm[jlm(l+1,mm-1)]
      Y.L.o[jlm(l,mm)]<-Ylm[jlm(l+1,mm  )]
      Y.L.M[jlm(l,mm)]<-Ylm[jlm(l+1,mm+1)]
#------------------------------------------------
      c.l[jlm(l,mm)]<-cl(l,mm)
      l.l[jlm(l,mm)]<-ll(l,mm)
#------------------------------------------------
      c.p[jlm(l,mm)]<-cp(l,mm)
      c.m[jlm(l,mm)]<-cm(l,mm)
      c.z[jlm(l,mm)]<-cz(l,mm)
#------------------------------------------------
      KmLm[jlm(l,mm)]<-Km(l+1,mm,mm-1)
      KoLo[jlm(l,mm)]<-Ko(l+1,mm,mm  )
      KpLM[jlm(l,mm)]<-Kp(l+1,mm,mm+1)
#------------------------------------------------
      KmlM[jlm(l,mm)]<-Km(l,mm,mm+1)
      Kolo[jlm(l,mm)]<-Ko(l,mm,mm  )
      Kplm[jlm(l,mm)]<-Kp(l,mm,mm-1)
#------------------------------------------------
      Km.lm[jlm(l,mm)]<-KYmlm(l,mm)
      Ko.lm[jlm(l,mm)]<-KYolm(l,mm)
      Kp.lm[jlm(l,mm)]<-KYplm(l,mm)
#------------------------------------------------
      Km.Lm[jlm(l,mm)]<-KYmLm(l,mm)
      Ko.Lm[jlm(l,mm)]<-KYoLm(l,mm)
      Kp.Lm[jlm(l,mm)]<-KYpLm(l,mm)
#------------------------------------------------
   }
#-------------------------------------------------------------------------------
# VECTOR SPHERICAL HARMONICS
#-------------------------------------------------------------------------------
   # R VERSOR (\bm{{\rm\hat{r}}})
   rm<-rep(sth*(cph-1i*sph)/sqrt(2),LMAX)
   rz<-rep(cth,LMAX)
   rp<-rep(sth*(cph+1i*sph)/sqrt(2),LMAX)
#------------------------------------------------
   # X - OK
   Xm<-l.l*c.m*Y.o.m/sqrt(2)
   Xz<-l.l*c.z*Y.o.o
   Xp<-l.l*c.p*Y.o.M/sqrt(2)
#------------------------------------------------
   # Y - OK
   Ym<-(-Kplm*Y.l.m+KmLm*Y.L.m)/sqrt(2)
   Yz<-  Kolo*Y.l.o+KoLo*Y.L.o
   Yp<-( KmlM*Y.l.M-KpLM*Y.L.M)/sqrt(2)
#------------------------------------------------
   # Y DEFINICAO
   Ym.o<-rm*Ylm[1:LMAX]
   Yz.o<-rz*Ylm[1:LMAX]
   Yp.o<-rp*Ylm[1:LMAX]
#------------------------------------------------
   # V
   Vm<-l.l*(-(c.l+1)*Kplm*Y.l.m-c.l*KmLm*Y.L.m)/sqrt(2)
   Vz<-l.l*( (c.l+1)*Kolo*Y.l.o-c.l*KoLo*Y.L.o)
   Vp<-l.l*( (c.l+1)*KmlM*Y.l.M+c.l*KpLM*Y.L.M)/sqrt(2)
#------------------------------------------------
   # V DEFINICAO
   Vm.o<-rm*Xz-rz*Xm
   Vz.o<-rp*Xm-rm*Xp
   Vp.o<-rz*Xp-rp*Xz
#------------------------------------------------
   # Yl,l-1,m 
   Ym.lm<--Km.lm*Y.l.m
   Yz.lm<- Ko.lm*Y.l.o
   Yp.lm<- Kp.lm*Y.l.M
#------------------------------------------------
   # Yl,l+1,m 
   Ym.Lm<--Km.Lm*Y.L.m
   Yz.Lm<--Ko.Lm*Y.L.o
   Yp.Lm<- Kp.Lm*Y.L.M 
#------------------------------------------------
   # Y - OK
   Ym.y<-(sqrt(c.l)*Ym.lm-sqrt(c.l+1)*Ym.Lm)/sqrt(2*c.l+1)
   Yz.y<-(sqrt(c.l)*Yz.lm-sqrt(c.l+1)*Yz.Lm)/sqrt(2*c.l+1)
   Yp.y<-(sqrt(c.l)*Yp.lm-sqrt(c.l+1)*Yp.Lm)/sqrt(2*c.l+1)
#------------------------------------------------
   # V - OK
   Vm.y<-(sqrt(c.l+1)*Ym.lm+sqrt(c.l)*Ym.Lm)/sqrt(2*c.l+1)
   Vz.y<-(sqrt(c.l+1)*Yz.lm+sqrt(c.l)*Yz.Lm)/sqrt(2*c.l+1)
   Vp.y<-(sqrt(c.l+1)*Yp.lm+sqrt(c.l)*Yp.Lm)/sqrt(2*c.l+1)
#-------------------------------------------------------------------------------
   # Y_{l,l-1}^m - OK
   Ym.lm.o<-(sqrt(c.l+1)*Vm+sqrt(c.l)*Ym)/sqrt(2*c.l+1) #OK
   Yz.lm.o<-(sqrt(c.l+1)*Vz+sqrt(c.l)*Yz)/sqrt(2*c.l+1) #OK
   Yp.lm.o<-(sqrt(c.l+1)*Vp+sqrt(c.l)*Yp)/sqrt(2*c.l+1) #OK
#-------------------------------------------------------------------------------
   # Y_{l,l+1}^m - OK
   Ym.Lm.o<-(sqrt(c.l)*Vm-sqrt(c.l+1)*Ym)/sqrt(2*c.l+1) #OK
   Yz.Lm.o<-(sqrt(c.l)*Vz-sqrt(c.l+1)*Yz)/sqrt(2*c.l+1) #OK
   Yp.Lm.o<-(sqrt(c.l)*Vp-sqrt(c.l+1)*Yp)/sqrt(2*c.l+1) #OK
#-------------------------------------------------------------------------------
# RESULTADOS
#-------------------------------------------------------------------------------
   Ylm<-Ylm[1:LMAX]
   Ym[1]<-Ylm[1]*sth*(cph-1i*sph)/sqrt(2)
   Yz[1]<-Ylm[1]*cth
   Yp[1]<-Ylm[1]*sth*(cph+1i*sph)/sqrt(2)
   U<-data.frame(l=c.l,m=c.z,
#                  r=l.l,cm=c.m,cp=c.p,
#                  KmLm,KoLo,KpLM,KmlM,Kolo,Kplm,
                  Ylm=Ylm,Xm,Xz,Xp,Ym,Yz,Yp,Vm,Vz,Vp,
                  Ym.lm,Yz.lm,Yp.lm,Ym.Lm,Yz.Lm,Yp.Lm,
                  Ym.y,Yz.y,Yp.y,Vm.y,Vz.y,Vp.y,
                  Ym.o,Yz.o,Yp.o,Vm.o,Vz.o,Vp.o,
                  Ym.lm.o,Yz.lm.o,Yp.lm.o,Ym.Lm.o,Yz.Lm.o,Yp.Lm.o)
   return(U)
}
#-------------------------------------------------------------------------------
# MULIPOLOS DE HANSEN
#-------------------------------------------------------------------------------
HansenMultipoles<-function(k,x,y,z,lmax){
   library(gsl)
#------------------------------------------------
   LMAX=lmax*(lmax+2)+1
   u<-SphericalHarmonics(x,y,z,lmax) 
   r<-sqrt(x^2+y^2+z^2)
   jl<-bessel_jl_steed_array(lmax+1,k*r) # Funcoes de bessel de 0 a lmax+1
#------------------------------------------------
   jl.m<-0 # Valor correto: cos(k*r)/(k*r); entretanto, M,N,L comecam em 1
   jl.o<-jl[1]
   jl.M<-jl[2]
#------------------------------------------------
   for(l in 1:lmax){
      n.rep<-2*l+1
      jl.m<-c(jl.m,rep(jl[l  ],n.rep))
      jl.o<-c(jl.o,rep(jl[l+1],n.rep))
      jl.M<-c(jl.M,rep(jl[l+2],n.rep))
   }
#------------------------------------------------
   # Constantes
   K<-2*u$l+1
   Kl<-u$l/K
   KL<-(u$l+1)/K
   Kll<-sqrt(u$l*(u$l+1))
   kl<-sqrt(Kl)
   kL<-sqrt(KL)
   p1<-KL*jl.m-Kl*jl.M # Derivada de jl(x)
   p2<-(jl.m+jl.M)/K   # jl(x)/x
#------------------------------------------------
   # M
   M.m<-jl.o*u$Xm
   M.z<-jl.o*u$Xz
   M.p<-jl.o*u$Xp
#------------------------------------------------
   # N
   N.m<-p1*u$Vm+Kll*p2*u$Ym
   N.z<-p1*u$Vz+Kll*p2*u$Yz
   N.p<-p1*u$Vp+Kll*p2*u$Yp
#------------------------------------------------
   # N_{lm} por Y_{j,l}^m
   Ny.m<-kL*jl.m*u$Ym.lm-kl*jl.M*u$Ym.Lm
   Ny.z<-kL*jl.m*u$Yz.lm-kl*jl.M*u$Yz.Lm
   Ny.p<-kL*jl.m*u$Yp.lm-kl*jl.M*u$Yp.Lm
#------------------------------------------------
   U<-data.frame(M.m,M.z,M.p,N.m,N.z,N.p)
#                 ,Ny.m,Ny.z,Ny.p)
}
