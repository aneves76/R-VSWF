#-------------------------------------------------------------------------------
# CALCULO EM LARGA ESCALA
#-------------------------------------------------------------------------------
source("RectangularWaveGuide.r")
source("PVectorSphericalWaveFunctions.r")
setEPS()
#-------------------------------------------------------------------------------
# Guia de onda cilindrico
RWG<-RectangularWaveGuide(TE=TRUE,kx,ky,kz,x+xo,y+yo,ze+zo)
#
tem.rwg<-array(RWG$Em,c(NP,NP,1))[,,1]
tez.rwg<-array(RWG$Ez,c(NP,NP,1))[,,1]
tep.rwg<-array(RWG$Ep,c(NP,NP,1))[,,1]
thm.rwg<-array(RWG$Hm,c(NP,NP,1))[,,1]
thz.rwg<-array(RWG$Hz,c(NP,NP,1))[,,1]
thp.rwg<-array(RWG$Hp,c(NP,NP,1))[,,1]
# Parte Real
pdf("ReEmRwgWfd.pdf");image(x+xo,y+yo,Re(tem.rwg),main="RWG Re Em",col=cm.colors(1024));grid();dev.off()
pdf("ReEzRwgWfd.pdf");image(x+xo,y+yo,Re(tez.rwg),main="RWG Re Ez",col=cm.colors(1024));grid();dev.off()
pdf("ReEpRwgWfd.pdf");image(x+xo,y+yo,Re(tep.rwg),main="RWG Re Ep",col=cm.colors(1024));grid();dev.off()
pdf("ReHmRwgWfd.pdf");image(x+xo,y+yo,Re(thm.rwg),main="RWG Re Hm",col=cm.colors(1024));grid();dev.off()
pdf("ReHzRwgWfd.pdf");image(x+xo,y+yo,Re(thz.rwg),main="RWG Re Hz",col=cm.colors(1024));grid();dev.off()
pdf("ReHpRwgWfd.pdf");image(x+xo,y+yo,Re(thp.rwg),main="RWG Re Hp",col=cm.colors(1024));grid();dev.off()
# Parte Imaginaria
pdf("ImEmRwgWfd.pdf");image(x+xo,y+yo,Im(tem.rwg),main="RWG Im Em",col=cm.colors(1024));grid();dev.off()
pdf("ImEzRwgWfd.pdf");image(x+xo,y+yo,Im(tez.rwg),main="RWG Im Ez",col=cm.colors(1024));grid();dev.off()
pdf("ImEpRwgWfd.pdf");image(x+xo,y+yo,Im(tep.rwg),main="RWG Im Ep",col=cm.colors(1024));grid();dev.off()
pdf("ImHmRwgWfd.pdf");image(x+xo,y+yo,Im(thm.rwg),main="RWG Im Hm",col=cm.colors(1024));grid();dev.off()
pdf("ImHzRwgWfd.pdf");image(x+xo,y+yo,Im(thz.rwg),main="RWG Im Hz",col=cm.colors(1024));grid();dev.off()
pdf("ImHpRwgWfd.pdf");image(x+xo,y+yo,Im(thp.rwg),main="RWG Im Hp",col=cm.colors(1024));grid();dev.off()
#-------------------------------------------------------------------------------
# Partial Wave Expansion
UVS<-PVectorSphericalWaveFunctions(k,x,y,ze,lmax,u$GTE,u$GTM)
#
tem.uvs<-array(UVS$Em,c(NP,NP,1))[,,1]
tez.uvs<-array(UVS$Ez,c(NP,NP,1))[,,1]
tep.uvs<-array(UVS$Ep,c(NP,NP,1))[,,1]
thm.uvs<-array(UVS$Hm,c(NP,NP,1))[,,1]
thz.uvs<-array(UVS$Hz,c(NP,NP,1))[,,1]
thp.uvs<-array(UVS$Hp,c(NP,NP,1))[,,1]
# Parte Real
pdf(paste("ReEmRwgPwe",lmax,".pdf",sep=""));image(x+xo,y+yo,Re(tem.uvs),main=paste("PWE Re Em",lmax),col=cm.colors(1024));grid();dev.off()
pdf(paste("ReEzRwgPwe",lmax,".pdf",sep=""));image(x+xo,y+yo,Re(tez.uvs),main=paste("PWE Re Ez",lmax),col=cm.colors(1024));grid();dev.off()
pdf(paste("ReEpRwgPwe",lmax,".pdf",sep=""));image(x+xo,y+yo,Re(tep.uvs),main=paste("PWE Re Ep",lmax),col=cm.colors(1024));grid();dev.off()
pdf(paste("ReHmRwgPwe",lmax,".pdf",sep=""));image(x+xo,y+yo,Re(thm.uvs),main=paste("PWE Re Hm",lmax),col=cm.colors(1024));grid();dev.off()
pdf(paste("ReHzRwgPwe",lmax,".pdf",sep=""));image(x+xo,y+yo,Re(thz.uvs),main=paste("PWE Re Hz",lmax),col=cm.colors(1024));grid();dev.off()
pdf(paste("ReHpRwgPwe",lmax,".pdf",sep=""));image(x+xo,y+yo,Re(thp.uvs),main=paste("PWE Re Hp",lmax),col=cm.colors(1024));grid();dev.off()
# Parte Imaginaria                ,sep=""  
pdf(paste("ReEmRwgPwe",lmax,".pdf",sep=""));image(x+xo,y+yo,Im(tem.uvs),main=paste("PWE Im Em",lmax),col=cm.colors(1024));grid();dev.off()
pdf(paste("ReEzRwgPwe",lmax,".pdf",sep=""));image(x+xo,y+yo,Im(tez.uvs),main=paste("PWE Im Ez",lmax),col=cm.colors(1024));grid();dev.off()
pdf(paste("ReEpRwgPwe",lmax,".pdf",sep=""));image(x+xo,y+yo,Im(tep.uvs),main=paste("PWE Im Ep",lmax),col=cm.colors(1024));grid();dev.off()
pdf(paste("ReHmRwgPwe",lmax,".pdf",sep=""));image(x+xo,y+yo,Im(thm.uvs),main=paste("PWE Im Hm",lmax),col=cm.colors(1024));grid();dev.off()
pdf(paste("ReHzRwgPwe",lmax,".pdf",sep=""));image(x+xo,y+yo,Im(thz.uvs),main=paste("PWE Im Hz",lmax),col=cm.colors(1024));grid();dev.off()
pdf(paste("ReHpRwgPwe",lmax,".pdf",sep=""));image(x+xo,y+yo,Im(thp.uvs),main=paste("PWE Im Hp",lmax),col=cm.colors(1024));grid();dev.off()
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
