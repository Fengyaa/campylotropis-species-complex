#!/usr/bin/env Rscript
# Base-R plotting; optional data.table speeds input. Raw GEMMA or processed full table.
a <- commandArgs(trailingOnly=TRUE)
if(length(a)<3) stop('Usage: Rscript plot_gemma.R FULL_RESULTS TRAIT OUTDIR [max_background_points=200000]')
input<-a[1];trait<-a[2];out<-a[3];maxbg<-if(length(a)>=4) as.integer(a[4]) else 200000L
if(is.na(maxbg)||maxbg<1) stop('max_background_points must be positive')
if(dir.exists(out)) stop('Output directory exists: use a new directory')
con<-if(grepl('\\.gz$',input)) gzfile(input,'rt') else file(input,'rt')
hdr<-strsplit(readLines(con,n=1), '[[:space:]]+')[[1]];close(con)
need<-intersect(c('rs','ps','p_wald','beta','se','valid_test'),hdr)
if(!all(c('rs','ps','p_wald') %in% need)) stop('Need full GEMMA table with rs, ps, p_wald; do not use signals-only file')
if(requireNamespace('data.table',quietly=TRUE)){
 # fread gzip support can depend on installed R.utils; use base fallback if unavailable.
 d<-tryCatch(as.data.frame(data.table::fread(input,select=need,na.strings=c('NA','NaN','.'))),error=function(e) NULL)
}else d<-NULL
if(is.null(d)){
 classes<-ifelse(hdr %in% need,NA,'NULL');con<-if(grepl('\\.gz$',input)) gzfile(input,'rt') else file(input,'rt')
 d<-read.table(con,header=TRUE,colClasses=classes,check.names=FALSE,comment.char='',na.strings=c('NA','NaN','.'));close(con)
}
M<-nrow(d)
if(M==0||anyDuplicated(d$rs)) stop('Empty input or duplicate SNP IDs')
if(any(!grepl('^HiC_scaffold[0-9]+:[0-9]+(:.*)?$',d$rs))) stop('Unexpected SNP ID format')
chr<-as.integer(sub('^HiC_scaffold([0-9]+):.*$','\\1',d$rs));pos<-as.numeric(sub('^HiC_scaffold[0-9]+:([0-9]+)(:.*)?$','\\1',d$rs))
if(any(pos!=d$ps)||!setequal(unique(chr),1:11)) stop('ID/position mismatch or incomplete 11-chromosome input')
p<-as.numeric(d$p_wald)
if(any(p<0|p>1,na.rm=TRUE)) stop('P outside [0,1]')
ok<-is.finite(p)
if('valid_test' %in% names(d))ok<-ok & !is.na(d$valid_test)&as.logical(d$valid_test)
if('beta' %in% names(d))ok<-ok & is.finite(d$beta)
if('se' %in% names(d))ok<-ok & is.finite(d$se)&d$se>0
if(!any(ok))stop('No valid tests')
limits<-sapply(1:11,function(c)max(pos[chr==c]));gap<-max(limits)*.015
starts<-c(0,head(cumsum(limits+gap),-1));mid<-starts+limits/2
ids<-d$rs[ok];chr<-chr[ok];pos<-pos[ok];p<-p[ok];rm(d);gc(verbose=FALSE)
# Clamp only zero P for visualization; retain count and record plotting floor.
floorp<-1e-300;y<--log10(pmax(p,floorp));x<-pos+starts[chr]
bonf<-.05/M;sugg<-1/M;sig<-p<=bonf
bg<-which(!sig);if(length(bg)>maxbg)bg<-bg[unique(round(seq(1,length(bg),length.out=maxbg)))]
ix<-c(bg,which(sig));ix<-ix[order(x[ix])]
# QQ positions from ALL valid tests; reduce rendered points only.
ord<-order(p);n<-length(p);ranks<-unique(c(seq_len(min(1000,n)),round(exp(seq(log(1),log(n),length.out=5000))),n))
qx<--log10((ranks-.5)/n);qy<--log10(pmax(p[ord[ranks]],floorp))
lo<--log10(qbeta(.975,ranks,n-ranks+1));hi<--log10(qbeta(.025,ranks,n-ranks+1));hi<-pmin(hi,300)
# Descriptive chi-square transformation only, not proof of valid LMM calibration.
lambda<-median(qchisq(p,df=1,lower.tail=FALSE))/qchisq(.5,df=1)
dir.create(out,recursive=TRUE)
manhattan<-function(){
 # At final size: 11 pt title, 10 pt axis labels, 9 pt ticks and legend.
 par(mar=c(4,4.5,2.5,1),mgp=c(2.7,.7,0),tcl=-.25,cex.main=1.1,cex.lab=1,cex.axis=.9);plot(x[ix]/1e6,y[ix],pch=20,cex=.45,col=ifelse(sig[ix],'#D55E00',ifelse(chr[ix]%%2,'#546E7A','#A8B6BC')),xaxt='n',xlab='Chromosome',ylab=expression(-log[10](P)),main=paste(trait,'GEMMA association'),ylim=c(0,max(y,-log10(bonf))*1.08))
 axis(1,at=mid/1e6,labels=1:11);abline(h=-log10(bonf),col='#B22222',lty=2);abline(h=-log10(sugg),col='#355C9A',lty=3)
 legend('topright',c('0.05 / all reported tests','1 / all reported tests'),col=c('#B22222','#355C9A'),lty=c(2,3),bty='n',cex=.9)
}
qq<-function(){
 par(mar=c(5,5,3,1));plot(qx,qy,type='n',xlab=expression(Expected~~-log[10](P)),ylab=expression(Observed~~-log[10](P)),main=paste(trait,'QQ plot'),ylim=c(0,max(qy,hi,na.rm=TRUE)*1.04))
 polygon(c(qx,rev(qx)),c(lo,rev(hi)),col='#DCE5ED',border=NA);abline(0,1,col='grey40',lty=2);points(qx,qy,pch=20,cex=.55,col='#275D7C')
 mtext(sprintf('Valid tests: %s | descriptive lambda: %.3f',format(n,big.mark=','),lambda),side=3,line=.25,cex=.75)
}
for(kind in c('manhattan','qq')){
 fn<-if(kind=='manhattan')manhattan else qq
 # Manhattan size: 200 x 90 mm; QQ retains its original 6 x 6 inch size.
 width_in<-if(kind=='manhattan')200/25.4 else 6
 height_in<-if(kind=='manhattan')90/25.4 else 6
 pointsize<-if(kind=='manhattan')10 else 12
 png(file.path(out,paste0(trait,'.',kind,'.png')),width=width_in,height=height_in,units='in',res=300,pointsize=pointsize);fn();dev.off()
 pdf(file.path(out,paste0(trait,'.',kind,'.pdf')),width=width_in,height=height_in,pointsize=pointsize,useDingbats=FALSE);fn();dev.off()
}
write.table(data.frame(trait=trait,reported_tests=M,valid_tests=n,bonferroni=bonf,suggestive=sugg,significant=sum(sig),zero_P=sum(p==0),plot_floor=floorp,rendered_manhattan=length(ix),rendered_QQ=length(ranks),lambda_descriptive=lambda),file.path(out,'plot_summary.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
writeLines(c('Input must be ALL reported tests, never a significant-only subset.','Thresholds are per-trait; no correction across eight traits.','All significant SNPs drawn; nonsignificant background deterministically downsampled only for rendering.','QQ ranks and lambda use all valid P values; QQ displayed points are subsampled.','QQ ribbon is pointwise 95% under independent uniform P values; genomic LD violates independence, so the ribbon is illustrative, not a calibrated simultaneous test.','Lambda is a descriptive chi-square(1) transformation; do not interpret as proof of confounding or use it to automatically correct LMM P values.','Zero P values plotted at 1e-300; positive P values below that floor are also capped for display.','Chromosome offsets use maximum reported marker positions, not assembly chromosome lengths.'),file.path(out,'plot_notes.txt'))
cat('Saved plots to',out,'\n')
