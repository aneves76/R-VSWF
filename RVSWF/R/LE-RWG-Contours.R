#-------------------------------------------------------------------------------
# CALCULO EM LARGA ESCALA
#-------------------------------------------------------------------------------
setEPS()
#-------------------------------------------------------------------------------
# Guia de onda cilindrico
# Parte Real
pdf("ReTemRwgWfd.pdf");contour(x+xo,y+yo,Re(tem.rwg),main="RWG Re Em");grid();dev.off()
pdf("ReTezRwgWfd.pdf");contour(x+xo,y+yo,Re(tez.rwg),main="RWG Re Ez");grid();dev.off()
pdf("ReTepRwgWfd.pdf");contour(x+xo,y+yo,Re(tep.rwg),main="RWG Re Ep");grid();dev.off()
pdf("ReThmRwgWfd.pdf");contour(x+xo,y+yo,Re(thm.rwg),main="RWG Re Hm");grid();dev.off()
pdf("ReThzRwgWfd.pdf");contour(x+xo,y+yo,Re(thz.rwg),main="RWG Re Hz");grid();dev.off()
pdf("ReThpRwgWfd.pdf");contour(x+xo,y+yo,Re(thp.rwg),main="RWG Re Hp");grid();dev.off()
# Parte Imaginaria
pdf("ImTemRwgWfd.pdf");contour(x+xo,y+yo,Im(tem.rwg),main="RWG Im Em");grid();dev.off()
pdf("ImTezRwgWfd.pdf");contour(x+xo,y+yo,Im(tez.rwg),main="RWG Im Ez");grid();dev.off()
pdf("ImTepRwgWfd.pdf");contour(x+xo,y+yo,Im(tep.rwg),main="RWG Im Ep");grid();dev.off()
pdf("ImThmRwgWfd.pdf");contour(x+xo,y+yo,Im(thm.rwg),main="RWG Im Hm");grid();dev.off()
pdf("ImThzRwgWfd.pdf");contour(x+xo,y+yo,Im(thz.rwg),main="RWG Im Hz");grid();dev.off()
pdf("ImThpRwgWfd.pdf");contour(x+xo,y+yo,Im(thp.rwg),main="RWG Im Hp");grid();dev.off()
#-------------------------------------------------------------------------------
# Partial Wave Expansion
# Parte Real
pdf("ReTemRwgPwe.pdf");contour(x+xo,y+yo,Re(tem.uvs),main=paste("PWE Re Em",lmax));grid();dev.off()
pdf("ReTezRwgPwe.pdf");contour(x+xo,y+yo,Re(tez.uvs),main=paste("PWE Re Ez",lmax));grid();dev.off()
pdf("ReTepRwgPwe.pdf");contour(x+xo,y+yo,Re(tep.uvs),main=paste("PWE Re Ep",lmax));grid();dev.off()
pdf("ReThmRwgPwe.pdf");contour(x+xo,y+yo,Re(thm.uvs),main=paste("PWE Re Hm",lmax));grid();dev.off()
pdf("ReThzRwgPwe.pdf");contour(x+xo,y+yo,Re(thz.uvs),main=paste("PWE Re Hz",lmax));grid();dev.off()
pdf("ReThpRwgPwe.pdf");contour(x+xo,y+yo,Re(thp.uvs),main=paste("PWE Re Hp",lmax));grid();dev.off()
# Parte Imaginaria
pdf("ImTemRwgPwe.pdf");contour(x+xo,y+yo,Im(tem.uvs),main=paste("PWE Im Em",lmax));grid();dev.off()
pdf("ImTezRwgPwe.pdf");contour(x+xo,y+yo,Im(tez.uvs),main=paste("PWE Im Ez",lmax));grid();dev.off()
pdf("ImTepRwgPwe.pdf");contour(x+xo,y+yo,Im(tep.uvs),main=paste("PWE Im Ep",lmax));grid();dev.off()
pdf("ImThmRwgPwe.pdf");contour(x+xo,y+yo,Im(thm.uvs),main=paste("PWE Im Hm",lmax));grid();dev.off()
pdf("ImThzRwgPwe.pdf");contour(x+xo,y+yo,Im(thz.uvs),main=paste("PWE Im Hz",lmax));grid();dev.off()
pdf("ImThpRwgPwe.pdf");contour(x+xo,y+yo,Im(thp.uvs),main=paste("PWE Im Hp",lmax));grid();dev.off()
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
