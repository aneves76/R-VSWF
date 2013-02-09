#-------------------------------------------------------------------------------
# CALCULO EM LARGA ESCALA
#-------------------------------------------------------------------------------
# Guia de onda cilindrico
# Parte Real
pdf("ReTemCwgWfd.pdf");contour(x+xo,y+yo,Re(tem.cwg),main="CWG Em");grid();dev.off()
pdf("ReTezCwgWfd.pdf");contour(x+xo,y+yo,Re(tez.cwg),main="CWG Ez");grid();dev.off()
pdf("ReTepCwgWfd.pdf");contour(x+xo,y+yo,Re(tep.cwg),main="CWG Ep");grid();dev.off()
pdf("ReThmCwgWfd.pdf");contour(x+xo,y+yo,Re(thm.cwg),main="CWG Hm");grid();dev.off()
pdf("ReThzCwgWfd.pdf");contour(x+xo,y+yo,Re(thz.cwg),main="CWG Hz");grid();dev.off()
pdf("ReThpCwgWfd.pdf");contour(x+xo,y+yo,Re(thp.cwg),main="CWG Hp");grid();dev.off()
# Parte Imaginaria
pdf("ImTemCwgWfd.pdf");contour(x+xo,y+yo,Im(tem.cwg),main="CWG Em");grid();dev.off()
pdf("ImTezCwgWfd.pdf");contour(x+xo,y+yo,Im(tez.cwg),main="CWG Ez");grid();dev.off()
pdf("ImTepCwgWfd.pdf");contour(x+xo,y+yo,Im(tep.cwg),main="CWG Ep");grid();dev.off()
pdf("ImThmCwgWfd.pdf");contour(x+xo,y+yo,Im(thm.cwg),main="CWG Hm");grid();dev.off()
pdf("ImThzCwgWfd.pdf");contour(x+xo,y+yo,Im(thz.cwg),main="CWG Hz");grid();dev.off()
pdf("ImThpCwgWfd.pdf");contour(x+xo,y+yo,Im(thp.cwg),main="CWG Hp");grid();dev.off()
#-------------------------------------------------------------------------------
# Partial Wave Expansion
# Parte Real
pdf("ReTemCwgPwe.pdf");contour(x+xo,y+yo,Re(tem.uvs),main=paste("PWE Em",lmax));grid();dev.off()
pdf("ReTezCwgPwe.pdf");contour(x+xo,y+yo,Re(tez.uvs),main=paste("PWE Ez",lmax));grid();dev.off()
pdf("ReTepCwgPwe.pdf");contour(x+xo,y+yo,Re(tep.uvs),main=paste("PWE Ep",lmax));grid();dev.off()
pdf("ReThmCwgPwe.pdf");contour(x+xo,y+yo,Re(thm.uvs),main=paste("PWE Hm",lmax));grid();dev.off()
pdf("ReThzCwgPwe.pdf");contour(x+xo,y+yo,Re(thz.uvs),main=paste("PWE Hz",lmax));grid();dev.off()
pdf("ReThpCwgPwe.pdf");contour(x+xo,y+yo,Re(thp.uvs),main=paste("PWE Hp",lmax));grid();dev.off()
# Parte Imaginaria
pdf("ImTemCwgPwe.pdf");contour(x+xo,y+yo,Im(tem.uvs),main=paste("PWE Em",lmax));grid();dev.off()
pdf("ImTezCwgPwe.pdf");contour(x+xo,y+yo,Im(tez.uvs),main=paste("PWE Ez",lmax));grid();dev.off()
pdf("ImTepCwgPwe.pdf");contour(x+xo,y+yo,Im(tep.uvs),main=paste("PWE Ep",lmax));grid();dev.off()
pdf("ImThmCwgPwe.pdf");contour(x+xo,y+yo,Im(thm.uvs),main=paste("PWE Hm",lmax));grid();dev.off()
pdf("ImThzCwgPwe.pdf");contour(x+xo,y+yo,Im(thz.uvs),main=paste("PWE Hz",lmax));grid();dev.off()
pdf("ImThpCwgPwe.pdf");contour(x+xo,y+yo,Im(thp.uvs),main=paste("PWE Hp",lmax));grid();dev.off()
#-------------------------------------------------------------------------------
