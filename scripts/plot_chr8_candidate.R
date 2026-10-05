# Paths are relative to this archive, independent of the caller's working directory.
.script <- sub('^--file=', '', commandArgs(FALSE)[grepl('^--file=', commandArgs(FALSE))])
root <- normalizePath(file.path(dirname(.script), '..'))
setwd(root)
# Run from the complex_analysis project root. Output: 180 x 75 mm.
o <- 'data/gwas/regional'
out <- 'results/candidate';dir.create(out,recursive=TRUE,showWarnings=FALSE)
aff <- read.delim('data/local_affinity/bin_weighting_sensitivity.tsv')
g <- subset(read.delim(file.path(o,'local_GWAS.tsv.gz')),scaffold=='HiC_scaffold8')
lead <- subset(read.delim(file.path(o,'lead_GWAS.tsv')),scaffold=='HiC_scaffold8')
gt <- subset(read.delim(file.path(o,'lead_genotypes_210.tsv')),scaffold=='HiC_scaffold8' & !is.na(GL))
variants <- read.delim('data/gwas/protein/verified_consequences.tsv')
variants <- subset(variants,consequence=='missense')
variant_rows <- match(variants$position,g$ps)
stopifnot(!anyNA(variant_rows),!anyDuplicated(g$ps),all(is.finite(g$minus_log10_P[variant_rows])))
variants$y <- g$minus_log10_P[variant_rows]
variants$label <- paste0(variants$reference_AA,variants$protein_position,variants$alternate_AA)
cols <- c(YS26='#2676A5',SM35='#CE6D2D',ML5='#8257A6')
draw <- function(){
 par(fig=c(0,.5,0,1),mar=c(4.3,4.5,.8,.8),mgp=c(2.9,.65,0),ps=9,cex.axis=.85,cex.lab=.9,tcl=-.25)
 wins <- c('left_100kb','candidate_100kb','MYB_gene','association_core','right_100kb')
 plot(NA,xlim=c(.5,5.5),ylim=c(-.055,.065),xaxt='n',xlab='',ylab='Distance to CYU - distance to CPO')
 abline(h=0,lty=2,col='gray50')
 axis(1,at=1:5,labels=c('Left\n100 kb','Candidate\n100 kb','MYB\ngene','Association\ncore','Right\n100 kb'),cex.axis=.8)
 for(i in 1:5) for(j in 1:3){
  v<-subset(aff,chrom=='Chr8' & window==wins[i] & population==names(cols)[j])$equal_bin_delta
  xx<-i+(j-2)*.19
  segments(xx,min(v),xx,max(v),col=cols[j],lwd=1.4)
  points(xx,median(v),pch=19,col=cols[j],cex=.8)
 }
 legend('topright',legend=names(cols),col=cols,pch=19,bty='n',cex=.8)
 par(fig=c(.5,1,0,1),new=TRUE,mar=c(4.3,4.2,.8,.8),mgp=c(2.5,.65,0))
 plot(g$ps/1e6,g$minus_log10_P,pch=20,cex=.4,col='#426B84',xlab='Position (Mb)',ylab=expression(-log[10](P)),ylim=c(0,max(g$minus_log10_P)*1.18))
 abline(h=-log10(1.44104885299717e-8),lty=2,col='#B33B36')
 peak<--log10(lead$p_wald)
 points(lead$ps/1e6,peak,pch=23,bg='#D55E00',col='#8C3900',cex=.9)
 text(lead$ps/1e6,peak+1.2,paste0('Lead SNP\n',format(lead$ps,big.mark=',',scientific=FALSE)),pos=2,cex=.78,col='#9B3C00')
 # Highlight verified missense SNPs at their observed GWAS coordinates.
 for(k in seq_len(nrow(variants))){
  vx<-variants$position[k]/1e6;vy<-variants$y[k]
  segments(vx-.19,vy+1.2,vx-.015,vy,col='#A13279',lwd=.8)
  points(vx,vy,pch=21,bg='#A13279',col='white',cex=1.05,lwd=.6)
  text(vx-.20,vy+1.2,paste0('MYB ',variants$label[k],'\n',format(variants$position[k],big.mark=',',scientific=FALSE)),pos=2,offset=.2,cex=.70,col='#A13279')
 }
 par(fig=c(.824,.986,.58,.89),new=TRUE,ps=6,mar=c(2.0,2.5,.3,.2),mgp=c(1.5,.4,0),cex.axis=.9,cex.lab=.9)
 tab<-table(factor(gt$GL,levels=0:1),factor(gt$GT,levels=c('0/0','0/1','1/1','./.')))
 bx<-barplot(tab,beside=TRUE,col=c('#0072B2','#D55E00'),border=NA,ylim=c(0,240),axes=FALSE,ylab='Individuals',cex.names=.85)
 axis(2,at=c(0,50,100,150),las=1)
 text(bx,tab+9,labels=ifelse(tab>0,tab,''),cex=.85)
 legend('topright',c('GL absent','GL present'),fill=c('#0072B2','#D55E00'),border=NA,bty='n',cex=.85,x.intersp=.5,y.intersp=.9)
}
pdf(file.path(out,'Chr8_affinity_GWAS_combined.pdf'),width=180/25.4,height=75/25.4,pointsize=9,useDingbats=FALSE);draw();dev.off()
png(file.path(out,'Chr8_affinity_GWAS_combined.png'),width=180,height=75,units='mm',res=300,pointsize=9);draw();dev.off()
