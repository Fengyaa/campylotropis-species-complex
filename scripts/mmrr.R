# Paths are relative to this archive, independent of the caller's working directory.
.script <- sub('^--file=', '', commandArgs(FALSE)[grepl('^--file=', commandArgs(FALSE))])
root <- normalizePath(file.path(dirname(.script), '..'))
setwd(root)
#!/usr/bin/env Rscript
base <- 'data/spatial'
out <- 'results/mmrr';dir.create(out,recursive=TRUE,showWarnings=FALSE)
source('vendor/PopGenReport/lgrMMRR.r')
y <- as.matrix(read.csv('results/mantel/fst_matrix_linearized.csv',row.names=1,check.names=FALSE))
meta <- read.csv('results/mantel/population_groups.csv')
stopifnot(identical(rownames(y),meta$population),identical(rownames(y),colnames(y)))
files <- c(Geographic='allpop_distance_matrix.txt',Environment='allpop_allenv_PCA_matrix.txt',
 Precipitation='allpop_precip_PCA_matrix.txt',Temperature='allpop_temp_PCA_matrix.txt',
 Solar='allpop_rsrad_PCA_matrix.txt',Wind='allpop_rwind_PCA_matrix.txt',Vapor='allpop_rvapr_PCA_matrix.txt')
x <- lapply(files,function(f){m<-as.matrix(read.table(file.path(base,f),header=TRUE,check.names=FALSE));dimnames(m)<-dimnames(y);m})
for(m in c(list(y),x))stopifnot(identical(dim(m),c(33L,33L)),all(is.finite(m)),isSymmetric(m,tol=1e-8))
groups <- list(Total=1:33,Cma=which(meta$species_test=='Cma'),Cyu=which(meta$species_test=='Cyu'),Cpo=which(meta$species_test=='Cpo'))
# Standardize unique off-diagonal pairs within each group, then mirror.
standardize <- function(m){v<-m[lower.tri(m)];z<-as.numeric(scale(v));stopifnot(all(is.finite(z)));m[,]<-0;m[lower.tri(m)]<-z;m<-m+t(m);m}
# Keep the official row-wise lower-triangle order for the independent OLS check.
unfold_check <- function(m) unlist(lapply(2:nrow(m),function(i)m[i,seq_len(i-1)]),use.names=FALSE)
run_group <- function(gi){
 g <- names(groups)[gi];ix<-groups[[gi]]
 yy<-standardize(y[ix,ix]);xx<-lapply(x,function(m)standardize(m[ix,ix]))
 jobs<-c(lapply(names(xx),function(v)list(type='Univariate',variable=v,terms=v)),
         lapply(names(xx)[-1],function(v)list(type='Geography-adjusted',variable=v,terms=c('Geographic',v))))
 rows<-list();models<-list()
 for(j in seq_along(jobs)){
  job<-jobs[[j]];seed<-20261002L+gi*100L+j
  cache<-file.path(out,sprintf('fit_%s_%02d.rds',g,j))
  if(file.exists(cache)){ans<-readRDS(cache)}else{
   set.seed(seed)
   ans<-lgrMMRR(gen.mat=yy,cost.mats=xx[job$terms],nperm=9999)
   saveRDS(ans,cache)
  }
  tab<-ans$mmrr.tab
  df<-as.data.frame(lapply(xx[job$terms],unfold_check));df$Y<-unfold_check(yy)
  fit<-lm(Y~.,data=df);sf<-summary(fit)
  for(term in job$terms){
   a<-tab[tab$layer==term,]
   stopifnot(nrow(a)==1,abs(a$coefficient-coef(fit)[term])<1e-10,
     abs(a$tstatistic-sf$coefficients[term,'t value'])<1e-8)
   rows[[length(rows)+1L]]<-data.frame(group=g,model=job$type,variable=job$variable,
     term=term,beta=a$coefficient,t=a$tstatistic,p=a$tpvalue,n_populations=length(ix),
     n_pairs=choose(length(ix),2),permutations=9999,seed=seed)
  }
  r2<-tab$r2[is.finite(tab$r2)];fp<-tab$Fpvalue[is.finite(tab$Fpvalue)]
  stopifnot(length(r2)==1,abs(r2-sf$r.squared)<1e-10)
  rho<-if(length(job$terms)==2)cor(df[[1]],df[[2]])else NA_real_
  models[[j]]<-data.frame(group=g,model=job$type,variable=job$variable,
    r2=r2,adjusted_r2=sf$adj.r.squared,F=unname(sf$fstatistic[1]),p_model=fp,
    predictor_correlation=rho,VIF=if(is.na(rho))1 else 1/(1-rho^2),seed=seed)
  message(g,' ',j,'/13: ',job$type,' ',job$variable,' complete')
 }
 list(coefficients=do.call(rbind,rows),models=do.call(rbind,models))
}
# Input fingerprint safeguards the cached permutation results.
inputs<-c('results/mantel/fst_matrix_linearized.csv',
 'results/mantel/population_groups.csv',file.path(base,files),
 'vendor/PopGenReport/lgrMMRR.r',file.path(root,'scripts/mmrr.R'))
hash<-tools::md5sum(inputs);mf<-file.path(out,'input_md5.rds')
if(file.exists(mf))stopifnot(identical(hash,readRDS(mf)))else saveRDS(hash,mf)
started<-Sys.time()
res<-parallel::mclapply(seq_along(groups),run_group,mc.cores=as.integer(Sys.getenv("CAMPY_CORES","1")),mc.set.seed=FALSE)
stopifnot(!any(vapply(res,inherits,logical(1),'try-error')))
a<-do.call(rbind,lapply(res,`[[`,'coefficients'));b<-do.call(rbind,lapply(res,`[[`,'models'))
a$p_BH_all_coefficients<-p.adjust(a$p,'BH')
a$stars<-ifelse(a$p<.001,'***',ifelse(a$p<.01,'**',ifelse(a$p<.05,'*','ns')))
write.csv(a,file.path(out,'MMRR_coefficients.csv'),row.names=FALSE)
write.csv(b,file.path(out,'MMRR_models.csv'),row.names=FALSE)
main<-a[a$term==a$variable,];main$p_BH_52<-p.adjust(main$p,'BH')
write.csv(main,file.path(out,'MMRR_table_long.csv'),row.names=FALSE)
for(metric in c('beta','p','p_BH_52')){
 wide<-data.frame(variable=names(files))
 for(model in c('Univariate','Geography-adjusted'))for(g in names(groups)){
  rr<-main[main$model==model & main$group==g,]
  wide[[paste(model,g,sep='_')]]<-rr[[metric]][match(wide$variable,rr$variable)]
 }
 write.csv(wide,file.path(out,paste0('table_',metric,'.csv')),row.names=FALSE,na='')
}
write.table(data.frame(file=names(hash),md5=unname(hash)),file.path(out,'input_md5.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
capture.output(sessionInfo(),file=file.path(out,'sessionInfo.txt'))
writeLines(paste('Elapsed seconds:',as.numeric(difftime(Sys.time(),started,units='secs'))),file.path(out,'timing.txt'))
print(main[,c('group','model','variable','beta','p')],row.names=FALSE)
