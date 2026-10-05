# Paths are relative to this archive, independent of the caller's working directory.
.script <- sub('^--file=', '', commandArgs(FALSE)[grepl('^--file=', commandArgs(FALSE))])
root <- normalizePath(file.path(dirname(.script), '..'))
setwd(root)
o <- 'results/hybrid';dir.create(o,recursive=TRUE,showWarnings=FALSE)
g <- as.matrix(read.delim('data/hybrid/GT.tsv.gz',row.names=1,check.names=FALSE)); storage.mode(g)<-'numeric'
dp <- as.matrix(read.delim('data/hybrid/DP.tsv.gz',row.names=1,check.names=FALSE));storage.mode(dp)<-'numeric'
m <- read.delim('data/hybrid/metadata.tsv',stringsAsFactors=FALSE)
f <- read.delim(gzfile('data/hybrid/marker_frequencies.tsv.gz'));stopifnot(identical(f$ID,rownames(g)),identical(m$sample,colnames(g)))
keep <- readLines('data/hybrid/pooled_ref_ld50kb.prune.in');main <- which(f$ID %in% keep)
p <- (f$p_ALT_CPO*40*f$call_CPO+.5)/(40*f$call_CPO+1);q <- (f$p_ALT_CYU*40*f$call_CYU+.5)/(40*f$call_CYU+1)
# Columns are lineage-pair probabilities: CPO/CPO, CPO/CYU, CYU/CYU.
prob <- function(v,p,q) {
 P<-cbind(ifelse(v==0,(1-p)^2,ifelse(v==1,2*p*(1-p),p^2)),ifelse(v==0,(1-p)*(1-q),ifelse(v==1,p*(1-q)+q*(1-p),p*q)),ifelse(v==0,(1-q)^2,ifelse(v==1,2*q*(1-q),q^2)))
 P
}
classes<-rbind(CPO=c(1,0,0),CYU=c(0,0,1),F1=c(0,1,0),F2=c(.25,.5,.25),BC_CPO=c(.5,.5,0),BC_CYU=c(0,.5,.5))
fit <- function(v,p,q) {
 good<-!is.na(v);P<-prob(v[good],p[good],q[good]);n<-nrow(P)
 if(n<10)return(c(S=NA,H=NA,n=n,raw_heter=NA,logLik=NA,F1_loss=NA))
 ll<-function(w)sum(log(pmax(as.vector(P%*%w),1e-300)))
 soft<-function(z){a<-c(z,0);a<-exp(a-max(a));a/sum(a)}
 ww<-function(z)c(z,1-sum(z))
 objective<-function(z)-ll(ww(z))
 grad<-function(z){den<-pmax(as.vector(P%*%ww(z)),1e-300);-colSums((P[,1:2,drop=FALSE]-P[,3])/den)}
 opt<-constrOptim(c(1/3,1/3),objective,grad,ui=rbind(c(1,0),c(0,1),c(-1,-1)),ci=c(0,0,-1),mu=1e-7,control=list(reltol=1e-10,maxit=300))
 ws<-rbind(ww(opt$par),diag(3))
 for(pair in list(c(1,2),c(1,3),c(2,3))) {
  wfun<-function(t){w<-numeric(3);w[pair]<-c(t,1-t);w}
  z<-optimize(function(t)-ll(wfun(t)),c(0,1),tol=1e-8);ws<-rbind(ws,wfun(z$minimum))
 }
 ls<-apply(ws,1,ll);w<-ws[which.max(ls),];c(S=w[1]+w[2]/2,H=w[2],n=n,raw_heter=mean(v[good]==1),logLik=max(ls),F1_loss=max(ls)-ll(classes['F1',]))
}
run <- function(ids,mat=g,pp=p,qq=q) t(vapply(seq_len(ncol(mat)),function(i)fit(mat[ids,i],pp[ids],qq[ids]),numeric(6)))
a<-run(main);res<-cbind(m,as.data.frame(a));res$call_rate<-res$n/length(main);res$interpretation_scope<-ifelse(res$K3_CMA>.1 | res$morph_group=='CMA','CMA_or_three_lineage_model_not_valid',ifelse(res$call_rate<.8,'low_call_rate','CPO_CYU_exploratory'))
write.table(res,file.path(o,'individual_estimates.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
sens<-list(all_candidates=seq_len(nrow(g)),delta08=which(f$delta>=.8 & seq_len(nrow(g))%in%main))
# Greedy physical spacing is an additional comparison, not proof of linkage independence.
spacing<-integer();last<-list()
for(i in seq_len(nrow(g))){c<-f$CHROM[i];if(is.null(last[[c]])||f$POS[i]-last[[c]]>=50000){spacing<-c(spacing,i);last[[c]]<-f$POS[i]}}
sens$spacing50kb<-spacing
sr<-list();for(nm in names(sens)){b<-run(sens[[nm]]);sr[[nm]]<-cbind(m[,1:2],panel=nm,as.data.frame(b))}
gcap<-g;caps<-apply(dp,2,median,na.rm=TRUE)*3;gcap[dp>matrix(caps,nrow(dp),ncol(dp),byrow=TRUE)]<-NA
sr$depth<-cbind(m[,1:2],panel='DP_cap_3x_candidate_median',as.data.frame(run(main,gcap)))
write.table(do.call(rbind,sr),file.path(o,'sensitivity_estimates.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
target<-which(m$population%in%c('ML5','YS26','SM35'));set.seed(20260927)
# Leave-one-reference-population-out frequency sensitivity; fixed ascertained marker set.
lr<-list();for(pop in unique(m$population[m$reference_group!='none'])){
 i1<-which(m$reference_group=='CPO' & m$population!=pop);i2<-which(m$reference_group=='CYU' & m$population!=pop)
 pp<-(rowSums(g[,i1,drop=FALSE],na.rm=TRUE)+.5)/(2*rowSums(!is.na(g[,i1,drop=FALSE]))+1);qq<-(rowSums(g[,i2,drop=FALSE],na.rm=TRUE)+.5)/(2*rowSums(!is.na(g[,i2,drop=FALSE]))+1)
 b<-t(vapply(target,function(i)fit(g[main,i],pp[main],qq[main]),numeric(6)))
 lr[[pop]]<-cbind(m[target,1:2],omitted_population=pop,as.data.frame(b))
}
write.table(do.call(rbind,lr),file.path(o,'reference_sensitivity.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
# Conditional parametric calibration under independent loci and estimated parental frequencies.
sim<-list();for(cl in rownames(classes)){
 b<-t(replicate(100,{
  anc<-sample(1:3,length(main),replace=TRUE,prob=classes[cl,]);pa<-ifelse(anc==3,q[main],p[main]);pb<-ifelse(anc==1,p[main],q[main]);v<-rbinom(length(main),1,pa)+rbinom(length(main),1,pb)
  fit(v,p[main],q[main])
 }))
 sim[[cl]]<-data.frame(class=cl,b)
}
sim<-do.call(rbind,sim);write.table(sim,file.path(o,'simulated_calibration.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
# 1-Mb block resampling, conditional on chosen markers and parental frequencies.
blocks<-split(main,paste(f$CHROM[main],(f$POS[main]-1)%/%1000000));br<-list()
for(i in target){b<-t(replicate(100,{ids<-unlist(blocks[sample(seq_along(blocks),length(blocks),replace=TRUE)],use.names=FALSE);fit(g[ids,i],p[ids],q[ids])}));br[[as.character(i)]]<-data.frame(sample=m$sample[i],population=m$population[i],S_lo=quantile(b[,'S'],.025),S_hi=quantile(b[,'S'],.975),H_lo=quantile(b[,'H'],.025),H_hi=quantile(b[,'H'],.975))}
write.table(do.call(rbind,br),file.path(o,'target_block_bootstrap.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
writeLines(c(paste('Candidates',nrow(g)),paste('LD-pruned',length(main)),paste('delta08 pruned',length(sens$delta08)),paste('Spacing50kb',length(spacing)),paste('1Mb blocks',length(blocks))),file.path(o,'panel_counts.txt'))
print(res[target,c('sample','population','S','H','n','raw_heter','F1_loss')]);print(readLines(file.path(o,'panel_counts.txt')))
