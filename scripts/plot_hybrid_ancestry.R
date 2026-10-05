# Paths are relative to this archive, independent of the caller's working directory.
.script <- sub('^--file=', '', commandArgs(FALSE)[grepl('^--file=', commandArgs(FALSE))])
root <- normalizePath(file.path(dirname(.script), '..'))
setwd(root)
o <- 'results/hybrid'
r <- read.delim(file.path(o,'individual_estimates.tsv'))
s <- read.delim(file.path(o,'simulated_calibration.tsv'))
cols <- c(ML5='#A4478D',YS26='#247BA0',SM35='#D58031')
classcols <- c(CPO='#3A966B',CYU='#4E70AD',F1='#AC5396',F2='#888888',BC_CPO='#8BAA49',BC_CYU='#CB9C41')
draw <- function() {
 par(mar=c(3.1,3.25,.35,.35),mgp=c(1.9,.55,0),tcl=-.22,las=1,lwd=.65,cex.axis=.9,family='Arial')
 plot(NA,xlim=c(-.035,1.035),ylim=c(-.025,1.075),xaxs='i',yaxs='i',xlab='Hybrid index S (CPO ancestry)',ylab='Interclass heterozygosity H',cex.lab=1)
 lines(c(0,.5,1),c(0,1,0),col='gray55',lty=2,lwd=.8)
 valid <- r$interpretation_scope=='CPO_CYU_exploratory'
 points(r$S[valid],r$H[valid],pch=16,col='#B8BDC480',cex=.65)
 for(cl in names(classcols)) {z<-s[s$class==cl,];points(z$S,z$H,pch=16,cex=.3,col=adjustcolor(classcols[cl],.25))}
 for(pop in names(cols)) {z<-r[r$population==pop,];points(z$S,z$H,pch=21,bg=cols[pop],col='white',lwd=.45,cex=.9)}
 legend('topright',legend=names(cols),pch=16,col=cols,bty='n',cex=.85,pt.cex=.85,inset=c(.005,.005),y.intersp=.95,x.intersp=.65)
 text(c(.06,.94),c(.055,.055),c('CYU','CPO'),cex=.8)
 text(.5,1.045,'F1',cex=.85)
 text(.5,.615,'F2',cex=.8)
 text(c(.19,.82),c(.59,.59),c('BC-CYU','BC-CPO'),cex=.75)
}
# 60 x 60 mm at final layout size; vector PDF plus 600 dpi PNG.
cairo_pdf(file.path(o,'hybrid_index_heterozygosity_left_60mm_Arial6pt.pdf'),width=60/25.4,height=60/25.4,pointsize=6,family='Arial');draw();dev.off()
png(file.path(o,'hybrid_index_heterozygosity_left_60mm_Arial6pt.png'),width=60,height=60,units='mm',res=600,pointsize=6,type='cairo',family='Arial');draw();dev.off()
