#-------------------------------------------------------------------------------
# CALCULO EM LARGA ESCALA
#-------------------------------------------------------------------------------
# Guia de onda cilindrico
#
pdf("ReTemBbzWfd.pdf");contour(x+xo,y+yo,Re(tem.BBZ),main="BBZ Em");grid();dev.off()
pdf("ReTezBbzWfd.pdf");contour(x+xo,y+yo,Re(tez.BBZ),main="BBZ Ez");grid();dev.off()
pdf("ReTepBbzWfd.pdf");contour(x+xo,y+yo,Re(tep.BBZ),main="BBZ Ep");grid();dev.off()
pdf("ReThmBbzWfd.pdf");contour(x+xo,y+yo,Re(thm.BBZ),main="BBZ Hm");grid();dev.off()
pdf("ReThzBbzWfd.pdf");contour(x+xo,y+yo,Re(thz.BBZ),main="BBZ Hz");grid();dev.off()
pdf("ReThpBbzWfd.pdf");contour(x+xo,y+yo,Re(thp.BBZ),main="BBZ Hp");grid();dev.off()
#
pdf("ImTemBbzWfd.pdf");contour(x+xo,y+yo,Im(tem.BBZ),main="BBZ Em");grid();dev.off()
pdf("ImTezBbzWfd.pdf");contour(x+xo,y+yo,Im(tez.BBZ),main="BBZ Ez");grid();dev.off()
pdf("ImTepBbzWfd.pdf");contour(x+xo,y+yo,Im(tep.BBZ),main="BBZ Ep");grid();dev.off()
pdf("ImThmBbzWfd.pdf");contour(x+xo,y+yo,Im(thm.BBZ),main="BBZ Hm");grid();dev.off()
pdf("ImThzBbzWfd.pdf");contour(x+xo,y+yo,Im(thz.BBZ),main="BBZ Hz");grid();dev.off()
pdf("ImThpBbzWfd.pdf");contour(x+xo,y+yo,Im(thp.BBZ),main="BBZ Hp");grid();dev.off()
#-------------------------------------------------------------------------------
# Partial Wave Expansion
# Parte Real
pdf("ReTemBbzPwe.pdf");contour(x+xo,y+yo,Re(tem.uvs),main=paste("PWE Em",lmax));grid();dev.off()
pdf("ReTezBbzPwe.pdf");contour(x+xo,y+yo,Re(tez.uvs),main=paste("PWE Ez",lmax));grid();dev.off()
pdf("ReTepBbzPwe.pdf");contour(x+xo,y+yo,Re(tep.uvs),main=paste("PWE Ep",lmax));grid();dev.off()
pdf("ReThmBbzPwe.pdf");contour(x+xo,y+yo,Re(thm.uvs),main=paste("PWE Hm",lmax));grid();dev.off()
pdf("ReThzBbzPwe.pdf");contour(x+xo,y+yo,Re(thz.uvs),main=paste("PWE Hz",lmax));grid();dev.off()
pdf("ReThpBbzPwe.pdf");contour(x+xo,y+yo,Re(thp.uvs),main=paste("PWE Hp",lmax));grid();dev.off()
# Parte Imaginaria
pdf("ImTemBbzPwe.pdf");contour(x+xo,y+yo,Im(tem.uvs),main=paste("PWE Em",lmax));grid();dev.off()
pdf("ImTezBbzPwe.pdf");contour(x+xo,y+yo,Im(tez.uvs),main=paste("PWE Ez",lmax));grid();dev.off()
pdf("ImTepBbzPwe.pdf");contour(x+xo,y+yo,Im(tep.uvs),main=paste("PWE Ep",lmax));grid();dev.off()
pdf("ImThmBbzPwe.pdf");contour(x+xo,y+yo,Im(thm.uvs),main=paste("PWE Hm",lmax));grid();dev.off()
pdf("ImThzBbzPwe.pdf");contour(x+xo,y+yo,Im(thz.uvs),main=paste("PWE Hz",lmax));grid();dev.off()
pdf("ImThpBbzPwe.pdf");contour(x+xo,y+yo,Im(thp.uvs),main=paste("PWE Hp",lmax));grid();dev.off()
#-------------------------------------------------------------------------------



