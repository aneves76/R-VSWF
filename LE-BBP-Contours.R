#-------------------------------------------------------------------------------
# CALCULO EM LARGA ESCALA
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Guia de onda cilindrico
# Parte Real
pdf("ReTemBbpWfd.pdf");contour(x+xo,y+yo,Re(tem.BBP),main="BBP Em");grid();dev.off()
pdf("ReTezBbpWfd.pdf");contour(x+xo,y+yo,Re(tez.BBP),main="BBP Ez");grid();dev.off()
pdf("ReTepBbpWfd.pdf");contour(x+xo,y+yo,Re(tep.BBP),main="BBP Ep");grid();dev.off()
pdf("ReThmBbpWfd.pdf");contour(x+xo,y+yo,Re(thm.BBP),main="BBP Hm");grid();dev.off()
pdf("ReThzBbpWfd.pdf");contour(x+xo,y+yo,Re(thz.BBP),main="BBP Hz");grid();dev.off()
pdf("ReThpBbpWfd.pdf");contour(x+xo,y+yo,Re(thp.BBP),main="BBP Hp");grid();dev.off()
# ParteTmaBnaWa
pdf("ImTemBbpWfd.pdf");contour(x+xo,y+yo,Im(tem.BBP),main="BBP Em");grid();dev.off()
pdf("ImTezBbpWfd.pdf");contour(x+xo,y+yo,Im(tez.BBP),main="BBP Ez");grid();dev.off()
pdf("ImTepBbpWfd.pdf");contour(x+xo,y+yo,Im(tep.BBP),main="BBP Ep");grid();dev.off()
pdf("ImThmBbpWfd.pdf");contour(x+xo,y+yo,Im(thm.BBP),main="BBP Hm");grid();dev.off()
pdf("ImThzBbpWfd.pdf");contour(x+xo,y+yo,Im(thz.BBP),main="BBP Hz");grid();dev.off()
pdf("ImThpBbpWfd.pdf");contour(x+xo,y+yo,Im(thp.BBP),main="BBP Hp");grid();dev.off()
#-------------------------------------------------------------------------------
# Partial Wave Expansion
# Parte Real
pdf("ReTemBbpPwe.pdf");contour(x+xo,y+yo,Re(tem.uvs),main=paste("PWE Em",lmax));grid();dev.off()
pdf("ReTezBbpPwe.pdf");contour(x+xo,y+yo,Re(tez.uvs),main=paste("PWE Ez",lmax));grid();dev.off()
pdf("ReTepBbpPwe.pdf");contour(x+xo,y+yo,Re(tep.uvs),main=paste("PWE Ep",lmax));grid();dev.off()
pdf("ReThmBbpPwe.pdf");contour(x+xo,y+yo,Re(thm.uvs),main=paste("PWE Hm",lmax));grid();dev.off()
pdf("ReThzBbpPwe.pdf");contour(x+xo,y+yo,Re(thz.uvs),main=paste("PWE Hz",lmax));grid();dev.off()
pdf("ReThpBbpPwe.pdf");contour(x+xo,y+yo,Re(thp.uvs),main=paste("PWE Hp",lmax));grid();dev.off()
# ParteTmaBnaWa
pdf("ImTemBbpPwe.pdf");contour(x+xo,y+yo,Im(tem.uvs),main=paste("PWE Em",lmax));grid();dev.off()
pdf("ImTezBbpPwe.pdf");contour(x+xo,y+yo,Im(tez.uvs),main=paste("PWE Ez",lmax));grid();dev.off()
pdf("ImTepBbpPwe.pdf");contour(x+xo,y+yo,Im(tep.uvs),main=paste("PWE Ep",lmax));grid();dev.off()
pdf("ImThmBbpPwe.pdf");contour(x+xo,y+yo,Im(thm.uvs),main=paste("PWE Hm",lmax));grid();dev.off()
pdf("ImThzBbpPwe.pdf");contour(x+xo,y+yo,Im(thz.uvs),main=paste("PWE Hz",lmax));grid();dev.off()
pdf("ImThpBbpPwe.pdf");contour(x+xo,y+yo,Im(thp.uvs),main=paste("PWE Hp",lmax));grid();dev.off()
#-------------------------------------------------------------------------------
