# Paths are relative to this archive, independent of the caller's working directory.
.script <- sub('^--file=', '', commandArgs(FALSE)[grepl('^--file=', commandArgs(FALSE))])
root <- normalizePath(file.path(dirname(.script), '..'))
setwd(root)
o <- 'results/local_affinity'
d <- read.delim(file.path(o,'bin_weighting_sensitivity.tsv'))
cols <- c(YS26='#2676A5',SM35='#CE6D2D',ML5='#8257A6')
draw <- function() {
par(mfrow=c(1,2),mar=c(6,5,3,1))
for(chr in c('Chr8','Chr6')) {
 wins <- if(chr=='Chr8') c('left_100kb','candidate_100kb','MYB_gene','association_core','right_100kb') else c('left_100kb','candidate_100kb','right_100kb')
 labs <- if(chr=='Chr8') c('Left\n100 kb','Candidate\n100 kb','MYB\ngene','Association\ncore','Right\n100 kb') else c('Left\n100 kb','Candidate\n100 kb','Right\n100 kb')
 plot(NA,xlim=c(.5,length(wins)+.5),ylim=c(-.055,.065),xaxt='n',xlab='',ylab='Distance to CYU - distance to CPO',main=paste(chr,'local genotype affinity'))
 abline(h=0,lty=2,col='gray50');axis(1,at=seq_along(wins),labels=labs,cex.axis=.8)
 for(i in seq_along(wins)) for(j in seq_along(cols)) {
  v <- subset(d,chrom==chr & window==wins[i] & population==names(cols)[j])$equal_bin_delta
  xx <- i+(j-2)*.19
  segments(xx,min(v),xx,max(v),col=cols[j],lwd=2)
  points(xx,median(v),pch=19,col=cols[j])
 }
 mtext('Positive: closer to CPO | Negative: closer to CYU',side=3,line=.2,cex=.7)
 legend('topright',legend=names(cols),col=cols,pch=19,bty='n',cex=.8)
}
mtext('Equal weight per occupied 1-kb bin. Points: median; bars: individual range (not confidence intervals).',side=1,outer=TRUE,line=-1.2,cex=.75)
}
pdf(file.path(o,'local_affinity.pdf'),width=12,height=5);draw();dev.off()
png(file.path(o,'local_affinity.png'),width=2400,height=1000,res=200);draw();dev.off()
