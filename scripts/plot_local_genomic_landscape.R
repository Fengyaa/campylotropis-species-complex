#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))
.script <- sub('^--file=', '', commandArgs(FALSE)[grepl('^--file=', commandArgs(FALSE))])
root <- normalizePath(file.path(dirname(.script), '..'));setwd(root)

here <- 'results/local';dir.create(here,recursive=TRUE,showWarnings=FALSE)
land <- 'data/windows'
w <- fread(file.path(land,'genomic_windows.tsv.gz'),na.strings=c('NA','nan','NaN',''))
w <- w[sites>=100]
setorder(w,chr,start)
target_chr <- 'HiC_scaffold11';lo <- 30000000;hi <- 50000000
loc <- w[chr==target_chr & end>lo & start<=hi]
stopifnot(nrow(loc)>0)
fields <- c('Fst_Cpo_Cyu','Mean_fst_sym','species_trio_fd_recomputed','Mean_f_d_sym_2026')
short <- c('Species FST','Sympatric FST','Species fd','Sympatric fd')
colors <- c('#315F78','#73A6B5','#BD673B','#E0A869')
sel <- list();records <- list()
for(j in seq_along(fields)){
 v <- w[[fields[j]]]
 eligible <- is.finite(v)
 if(j>=3)eligible <- eligible & v>0 & v<1
 for(p in c(1,2.5,5)){
  cut <- as.numeric(quantile(v[eligible],1-p/100,type=7))
  mark <- loc[is.finite(get(fields[j])) & get(fields[j])>=cut]
  records[[paste(j,p)]] <- data.table(statistic=fields[j],top_percent=p,genomewide_cutoff=cut,
    genomewide_eligible=sum(eligible),local_windows=nrow(mark),local_100kb_bins=uniqueN((mark$start-1)%/%100000))
  if(p==1){mark[,`:=`(statistic=fields[j],top_percent=p)];sel[[j]] <- mark}
 }
}
fwrite(rbindlist(records),file.path(here,'threshold_comparison.tsv'),sep='\t')
marks <- rbindlist(sel,fill=TRUE)
fwrite(marks[,.(chr,start,end,statistic,top_percent)],file.path(here,'top1_candidate_windows.tsv'),sep='\t')
g <- fread('data/annotations/gene_coordinates.tsv.gz',header=FALSE,col.names=c('chr','gene','start','end'))[chr==target_chr & start<=hi & end>=lo]
nbs <- unique(sub('[.]m[0-9]+$','',trimws(readLines('data/annotations/NBS_gene_ids.txt'))))
ng <- g[gene%in%nbs]
stopifnot(nrow(ng)==10)
fwrite(ng,file.path(here,'local_NBS_genes.tsv'),sep='\t')
rp <- fread('data/annotations/repeat_union_100kb_chr11.tsv.gz')[chr==target_chr & end>lo & start<=hi]
stopifnot(nrow(rp)==200,all(rp$repeat_content>=0 & rp$repeat_content<=1))
# Variant-only pi and Dxy are retained as contextual SNP-set measures.
tracks <- list(
 list(label=expression(F[ST]~'species'),field='Fst_Cpo_Cyu',max=.30,color=colors[1],kind='line'),
 list(label=expression(F[ST]~'sympatric'),field='Mean_fst_sym',max=.30,color=colors[2],kind='line'),
 list(label=expression(f[d]~'species'),field='species_trio_fd_recomputed',max=1,color=colors[3],kind='line'),
 list(label=expression(f[d]~'sympatric'),field='Mean_f_d_sym_2026',max=1,color=colors[4],kind='line'),
 list(label=expression(D[XY]~'species'),field='dxy_Cpo_Cyu',max=.25,color='#556278',kind='line'),
 list(label=expression(D[XY]~'sympatric'),field='Mean_dxy_sym',max=.25,color='#939DB0',kind='line'),
 list(label=expression(pi~'Cpo / Cyu'),field=NULL,max=.22,color=NA,kind='pi'),
 list(label='Gene positions',field=NULL,max=1,color=NA,kind='genes'),
 list(label='Repeat fraction',field=NULL,max=1,color=NA,kind='repeat'))

draw <- function(){
 par(mar=c(0,0,0,0),family='Helvetica',xpd=NA)
 plot.new();plot.window(xlim=c(27.35,50.8),ylim=c(.55,13.6),xaxs='i',yaxs='i')
 text(27.5,13.30,'Chr11: 30-50 Mb | local differentiation and introgression',adj=0,font=2,cex=.86)
 text(27.5,13.02,'Candidate marks: genome-wide top 1% of each statistic',adj=0,cex=.60,col='#586773')
 xs <- c(30,35,40,45,50)
 # Four separate candidate bands, with a minimum visible tick width.
 for(j in 1:4){
  top <- 12.71-(j-1)*.205;bot <- top-.15
  rect(30,bot,50,top,col='#F5F6F7',border=NA)
  a <- marks[statistic==fields[j]]
  x1 <- pmax(lo,a$start)/1e6;x2 <- pmin(hi,a$end)/1e6
  xm <- (x1+x2)/2
  rect(pmax(30,xm-.018),bot,pmin(50,xm+.018),top,col=colors[j],border=NA)
  text(29.77,(top+bot)/2,short[j],adj=1,col=colors[j],cex=.52)
 }
 for(i in seq_along(tracks)){
  t <- tracks[[i]];top <- 11.56-(i-1)*1.075;bottom <- top-.83
  scale <- function(v) bottom+pmax(0,pmin(v,t$max))/t$max*(top-bottom)
  rect(30,bottom,50,top,col=if(i%%2==1)'#F8F9FA' else '#F3F6F8',border=NA)
  segments(xs,bottom,xs,top,col='#E7EAEC',lwd=.45)
  # Faint NBS guides connect annotated positions across all local tracks.
  segments((ng$start+ng$end)/2e6,bottom,(ng$start+ng$end)/2e6,top,
           col=adjustcolor('#673C91',alpha.f=.20),lwd=.6)
  segments(30,bottom,50,bottom,col='#89939B',lwd=.42)
  text(27.5,(top+bottom)/2,t$label,adj=0,cex=.64,col=if(i<=4)t$color else '#41515E')
  if(t$kind=='line'){
   a <- loc[,.(x=(start+end)/2e6,y=get(t$field))]
   if(i<=4)a[,y:=pmax(y,0)]
   lines(a$x,scale(a$y),col=t$color,lwd=.42)
   points(a$x,scale(a$y),pch=16,cex=.075,col=t$color)
  }
  if(t$kind=='pi'){
   for(k in 1:2){v <- loc[[c('pi_Cpo','pi_Cyu')[k]]]
    lines((loc$start+loc$end)/2e6,scale(v),col=c('#B8743D','#7569AA')[k],lwd=.44)
   }
   text(49.7,top-.09,'Cpo',adj=1,col='#B8743D',cex=.42)
   text(49.7,top-.27,'Cyu',adj=1,col='#7569AA',cex=.42)
  }
  if(t$kind=='genes'){
   rect(pmax(lo,g$start)/1e6,bottom+.18,pmin(hi,g$end)/1e6,bottom+.29,col='#A7AFB5',border=NA)
   rect(pmax(lo,ng$start)/1e6,bottom+.14,pmin(hi,ng$end)/1e6,bottom+.48,col='#673C91',border=NA)
   # Genomic clusters identified by the existing BLAST candidate list.
   text(c(35.6,46.8),top-.12,c('NBS (4)','NBS (6)'),col='#673C91',cex=.49)
  }
  if(t$kind=='repeat'){
   rect(pmax(lo,rp$start)/1e6,bottom,pmin(hi,rp$end)/1e6,scale(rp$repeat_content),
        col='#AAA5B6',border=NA)
  }
  if(t$kind%in%c('line','pi','repeat')){
   text(29.75,c(bottom,top),c('0',format(t$max,trim=TRUE)),adj=1,cex=.42,col='#68737B')
  }
 }
 # Shared x axis with no duplicated per-track tick labels.
 segments(30,1.98,50,1.98,col='#738089',lwd=.6)
 segments(xs,1.98,xs,1.91,col='#738089',lwd=.6)
 text(xs,1.78,xs,cex=.62,col='#46515A')
 text(40,1.40,'Chromosome position (Mb)',cex=.72)
 text(27.5,.90,'Pi and Dxy are calculated on filtered SNP sites; repeat proportions use 100-kb bins.',
      adj=0,cex=.48,col='#586773')
}
for(fmt in c('pdf','png')){
 out <- file.path(here,paste0('Chr11_30_50Mb_local_top1_NBS_guides.',fmt))
 if(fmt=='pdf')pdf(out,width=190/25.4,height=185/25.4,useDingbats=FALSE)
 else png(out,width=190,height=185,units='mm',res=260,type='cairo')
 draw();dev.off()
}
inputs <- c(file.path(land,'genomic_windows.tsv.gz'),'data/annotations/repeat_union_100kb_chr11.tsv.gz',
            'data/annotations/gene_coordinates.tsv.gz','data/annotations/NBS_gene_ids.txt')
fwrite(data.table(path=inputs,md5=unname(tools::md5sum(inputs))),file.path(here,'input_manifest.tsv'),sep='\t')
cat('Saved Chr11 local plot: ',nrow(loc),'statistics windows, ',nrow(g),'genes, ',nrow(ng),'NBS.\n')
