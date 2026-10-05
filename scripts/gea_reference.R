# Paths are relative to this archive, independent of the caller's working directory.
.script <- sub('^--file=', '', commandArgs(FALSE)[grepl('^--file=', commandArgs(FALSE))])
root <- normalizePath(file.path(dirname(.script), '..'))
setwd(root)
suppressPackageStartupMessages(library(data.table))
out <- 'results/gea_reference';dir.create(out,recursive=TRUE,showWarnings=FALSE)
w <- fread('data/windows/genomic_windows.tsv.gz')
lens <- fread('data/annotations/chromosome_lengths.tsv')
setorder(w,chr,start);stopifnot(!anyDuplicated(w[,.(chr,start)]))
w[,id:=paste(chr,start)]
traits <- c('fd_MY1_CYU_CPO','fd_CY32_CPO_CYU','species_fd_common_support')
labels <- 'Z3.5_r0.5'
files <- setNames(file.path('data/gea',paste0(labels,'.shared_10kb_windows.tsv')),labels)
B <- 9999L;set.seed(20261003)
shifts <- lapply(lens$reference_length,function(n)sample.int(ceiling(n/10000),B,replace=TRUE));names(shifts)<-lens$chr
cross <- function(a,y) {
 n<-length(a);v<-Re(fft(Conj(fft(a))*fft(y),inverse=TRUE))/n
 stopifnot(abs(v[1]-sum(a*y))<1e-6*max(1,abs(sum(a*y))))
 for(s in c(1L,n%/%2L))stopifnot(abs(v[s+1]-sum(a*y[((seq_len(n)-1+s)%%n)+1]))<1e-6*max(1,abs(sum(a*y))))
 v
}
results<-list();qc<-list();cires<-list()
for(label in names(files)) {
 cand<-fread(files[[label]]);ids<-paste(cand$chromosome,cand$start);sel<-w$id%in%ids
 qc[[label]]<-data.table(candidate_set=label,input_windows=nrow(cand),present_in_parameter_table=sum(sel),sites_ge100=sum(sel&w$sites>=100))
 fwrite(w[sel & is.finite(sites) & sites>=100],file.path(out,paste0(label,'.candidate_parameters.tsv')),sep='\t')
 for(trait in traits) {
  y<-w[[trait]];ok<-is.finite(y)&is.finite(w$sites)&w$sites>=100
  if(trait%in%traits)ok<-ok&y>0&y<1
  if(grepl('^xpclr',trait)||trait=='Rho')ok<-ok&y>=0
  a<-sel&ok;b<-!sel&ok;stopifnot(sum(a)>1,sum(b)>1)
  obs<-mean(y[a])-mean(y[b]);ns<-nc<-numeric(B)
  for(cc in lens$chr){
   n<-ceiling(lens[chr==cc,reference_length]/10000);ix<-which(w$chr==cc);pos<-(w$start[ix]-1)%/%10000+1
   av<-ev<-yv<-numeric(n);av[pos]<-as.numeric(sel[ix]);ev[pos]<-as.numeric(ok[ix]);yv[pos[ok[ix]]]<-y[ix[ok[ix]]]
   ns<-ns+cross(av,yv)[shifts[[cc]]];nc<-nc+pmax(0,round(cross(av,ev)[shifts[[cc]]]))
  }
  stopifnot(all(nc>0&nc<sum(ok)))
  null<-ns/nc-(sum(y[ok])-ns)/(sum(ok)-nc)
  p<-min(1,2*min((1+sum(null>=obs))/(B+1),(1+sum(null<=obs))/(B+1)))
  results[[paste(label,trait)]]<-data.table(candidate_set=label,trait=trait,n_candidate=sum(a),n_background=sum(b),mean_candidate=mean(y[a]),mean_background=mean(y[b]),median_candidate=median(y[a]),median_background=median(y[b]),mean_difference=obs,relative_difference_percent=100*obs/mean(y[b]),standardized_difference=obs/sd(y[ok]),p_spatial=p,null_median=median(null))
  for(bp in c(1e6,5e6)){
   agg<-data.table(chr=w$chr,block=(w$start-1)%/%bp,sa=ifelse(a,y,0),na=as.integer(a),sb=ifelse(b,y,0),nb=as.integer(b))[,.(sa=sum(sa),na=sum(na),sb=sum(sb),nb=sum(nb)),by=.(chr,block)]
   totals<-matrix(0,1000,4)
   for(cc in unique(agg$chr)){
    mat<-as.matrix(agg[chr==cc,.(sa,na,sb,nb)]);k<-nrow(mat)
    weights<-t(rmultinom(1000,k,rep(1/k,k)));totals<-totals+weights%*%mat
   }
   dif<-totals[,1]/totals[,2]-totals[,3]/totals[,4];stopifnot(all(is.finite(dif)))
   ci<-quantile(dif,c(.025,.975))
   cires[[paste(label,trait,bp)]]<-data.table(candidate_set=label,trait=trait,block_bp=bp,difference_ci_low=ci[1],difference_ci_high=ci[2],standardized_ci_low=ci[1]/sd(y[ok]),standardized_ci_high=ci[2]/sd(y[ok]))
  }
 }
 cat('Completed',label,'\n')
}
r<-rbindlist(results);r[,q_BH_within_set:=p.adjust(p_spatial,'BH'),by=candidate_set];r[,q_BH_all_sets:=p.adjust(p_spatial,'BH')]
fwrite(r,file.path(out,'comparison_statistics.tsv'),sep='\t');fwrite(rbindlist(qc),file.path(out,'window_qc.tsv'),sep='\t');ci<-rbindlist(cires);fwrite(ci,file.path(out,'block_bootstrap_intervals.tsv'),sep='\t')
plotdata<-merge(r,ci[block_bp==1e6],by=c('candidate_set','trait'))
plotfun<-function(){
 par(mfrow=c(1,2),mar=c(4,12,3,1))
 for(label in labels){
  d<-plotdata[candidate_set==label];d<-d[match(traits,trait)];yy<-rev(seq_len(nrow(d)));xr<-range(c(d$standardized_ci_low,d$standardized_ci_high,0))
  plot(d$standardized_difference,yy,xlim=xr,ylim=c(.5,nrow(d)+.5),yaxt='n',ylab='',xlab='Mean difference / genome-wide SD',main=paste(label,'and BayPass'),pch=19,col=ifelse(d$q_BH_within_set<.05,'#B34D32','#647989'))
  abline(v=0,lty=2,col='gray60');segments(d$standardized_ci_low,yy,d$standardized_ci_high,yy,col='#647989')
  axis(2,at=yy,labels=d$trait,las=2,cex.axis=.7)
 }
}
pdf(file.path(out,'genetic_parameter_comparison.pdf'),width=14,height=9);plotfun();dev.off()
png(file.path(out,'genetic_parameter_comparison.png'),width=2800,height=1800,res=200);plotfun();dev.off()
fwrite(data.table(path=c(files,'data/windows/genomic_windows.tsv.gz'),md5=unname(tools::md5sum(c(files,'data/windows/genomic_windows.tsv.gz')))),file.path(out,'input_manifest.tsv'),sep='\t')
writeLines(capture.output(sessionInfo()),file.path(out,'sessionInfo.txt'))
print(r[candidate_set=='Z3.5_r0.5',.(trait,relative_difference_percent,p_spatial,q_BH_within_set)])
