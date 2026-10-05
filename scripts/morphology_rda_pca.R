#!/usr/bin/env Rscript
.script<-sub('^--file=','',commandArgs(FALSE)[grepl('^--file=',commandArgs(FALSE))])
root<-normalizePath(file.path(dirname(.script),'..'));setwd(root)
suppressPackageStartupMessages(library(vegan))
out<-'results/traits';dir.create(out,recursive=TRUE,showWarnings=FALSE)
write_tsv<-function(x,n)write.table(x,file.path(out,paste0(n,'.tsv')),sep='\t',quote=FALSE,row.names=FALSE)
d<-read.delim('data/traits/RDA_input_192.tsv',check.names=FALSE)
env<-read.delim('data/climate/population_climate_scaled.tsv',check.names=FALSE)
tt<-c('TLL','TLW','RL','PL','SL_SW','TLL_TLW','ST','GL')
cc<-c('bio3','bio7','bio10','bio13','bio15','srad_07')
stopifnot(nrow(d)==192,length(unique(d$code))==28,!anyDuplicated(d$sample))
for(v in cc){
 expected<-env[[if(v=='srad_07')'wc2.1_2.5m_srad_07' else v]][match(d$code,env$code)]
 stopifnot(all(is.finite(expected)),max(abs(expected-d[[v]]))<1e-10)
}
Y<-scale(as.matrix(d[,tt]));X<-as.data.frame(scale(d[,cc]))
m<-rda(Y~bio3+bio7+bio10+bio13+bio15+srad_07,data=X)
rs<-RsquareAdj(m);ev<-100*m$CCA$eig[1:2]/m$tot.chi
write_tsv(data.frame(n=nrow(d),populations=28,R2=rs$r.squared,adjusted_R2=rs$adj.r.squared,RDA1_percent_total=ev[1],RDA2_percent_total=ev[2]),'RDA_summary')
for(type in c('sites','bp','species')){
 z<-scores(m,display=type,scaling=2,choices=1:2)
 write_tsv(data.frame(id=if(type=='sites')d$sample else rownames(z),z),paste0('RDA_',type))
}
pop<-aggregate(d[,tt],list(code=d$code),mean)
ec<-grep('^bio[0-9]+$|^wc2[.]1_2[.]5m_(srad|wind|vapr)_[0-9]+$',names(env),value=TRUE)
stopifnot(length(ec)==55)
inputs<-list(population_traits=pop,individual_traits=d,environment=env)
vars<-list(population_traits=tt,individual_traits=tt,environment=ec);summary<-list()
for(nm in names(inputs)){
 x<-inputs[[nm]];fit<-prcomp(x[,vars[[nm]]],center=TRUE,scale.=TRUE);pc<-100*fit$sdev^2/sum(fit$sdev^2)
 summary[[nm]]<-data.frame(analysis=nm,n=nrow(x),PC1=pc[1],PC2=pc[2])
 write_tsv(data.frame(id=if('sample'%in%names(x))x$sample else x$code,fit$x),paste0(nm,'_PCA_scores'))
 write_tsv(data.frame(variable=rownames(fit$rotation),fit$rotation),paste0(nm,'_PCA_loadings'))
}
write_tsv(do.call(rbind,summary),'PCA_summary');print(rs);print(ev)
