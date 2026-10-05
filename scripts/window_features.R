# Paths are relative to this archive, independent of the caller's working directory.
.script <- sub('^--file=', '', commandArgs(FALSE)[grepl('^--file=', commandArgs(FALSE))])
root <- normalizePath(file.path(dirname(.script), '..'))
setwd(root)
here <- 'results/windows';dir.create(here,recursive=TRUE,showWarnings=FALSE)
suppressPackageStartupMessages({library(data.table);library(clusterProfiler);library(ggplot2)})
set.seed(20260927)
write_tsv <- function(x,n)fwrite(x,file.path(here,n),sep='\t',quote=FALSE,na='NA')
land <- file.path(root,'data/windows')
prev <- file.path(root,'data/annotations')
src <- file.path(land,'genomic_windows.tsv.gz')
w <- fread(src,na.strings=c('NA','nan','NaN',''))
annfile <- file.path(root,'data/annotations/window_union_coverage.tsv.gz')
ann <- fread(annfile)[,.(chr=chrom,start=start_0based+1,end=end_exclusive,gene_bp=gene_union_bp,repeat_bp=repeat_union_bp)]
stopifnot(!anyDuplicated(ann[,.(chr,start,end)]),!anyDuplicated(w[,.(chr,start,end)]))
w <- merge(w,ann,by=c('chr','start','end'),all.x=TRUE,sort=FALSE)
lens <- fread(file.path(root,'data/annotations/chromosome_lengths.tsv'))
w <- merge(w,lens,by='chr',all.x=TRUE,sort=FALSE)
w[,annotation_width:=pmin(end,reference_length)-start+1]
stopifnot(!anyNA(w$gene_bp),!anyNA(w$repeat_bp),all(w$annotation_width>0))
w[,`:=`(gene_coverage=gene_bp/annotation_width,repeat_coverage=repeat_bp/annotation_width)]
stopifnot(all(w$gene_coverage>=0&w$gene_coverage<=1),all(w$repeat_coverage>=0&w$repeat_coverage<=1))
w[,window_id:=paste(chr,start,end,sep=':')]
setorder(w,chr,start)
genes <- fread(file.path(root,'data/annotations/gene_coordinates.tsv.gz'),header=FALSE,col.names=c('chr','gene_id','start','end'))
setkey(genes,chr,start,end)
hits <- unique(foverlaps(w[,.(chr,start,end,window_id)],genes,type='any',nomatch=0L)[,.(window_id,gene_id)])
anno <- list(GO=fread(file.path(prev,'gene_GO.tsv.gz')),KEGG=fread(file.path(prev,'gene_KEGG.tsv.gz')))
contrasts <- c(FST_species='Fst_Cpo_Cyu',fd_species='species_trio_fd_recomputed',
 FST_sympatric='Mean_fst_sym',fd_sympatric='Mean_f_d_sym_2026')
labels <- c(FST_species='High FST (Cpo-Cyu)',fd_species='High fd (species trio)')
colors <- c(FST_species='#315F78',fd_species='#BD673B')
ps <- c(10,5,2.5,1)
traits <- c('gene_coverage','repeat_coverage','NBS_density')
nbsfile <- file.path(root,'data/annotations/NBS_gene_ids.txt')
nbs <- unique(sub('[.]m[0-9]+$','',trimws(readLines(nbsfile))))
nh <- hits[gene_id%in%nbs,.(NBS_count=uniqueN(gene_id)),by=window_id]
w <- merge(w,nh,by='window_id',all.x=TRUE,sort=FALSE)
w[is.na(NBS_count),NBS_count:=0L]
w[,NBS_density:=NBS_count/annotation_width*1e6]
setorder(w,chr,start)
trait_labels <- c('Gene\ncoverage','Repeat\ncoverage','Diversity\nC. polyantha','Diversity\nC. yunnanensis','Absolute divergence\n(DXY)')
thresholds <- list();selected_sets <- list();valid_sets <- list();enrichments <- list()

enrich <- function(selected,universe,db) {
  a <- anno[[db]];u <- intersect(universe,unique(a$gene_id));s <- intersect(selected,u)
  t2g <- unique(a[gene_id%in%u,.(ID,gene_id)])
  sets <- split(t2g$gene_id,t2g$ID);sets <- sets[lengths(sets)>=10&lengths(sets)<=500]
  h <- lapply(sets,intersect,s);k <- lengths(h);M <- lengths(sets);N <- length(u);n <- length(s)
  r <- data.table(ID=names(sets),Count=as.integer(k),BgCount=as.integer(M),SelectedAnnotatedGenes=n,BackgroundAnnotatedGenes=N,
    FoldEnrichment=(k/n)/(M/N),pvalue=phyper(k-1,M,N-M,n,lower.tail=FALSE),geneID=vapply(h,paste,collapse='/',FUN.VALUE=character(1)))
  r[,FDR_fixed_family:=p.adjust(pvalue,'BH')]
  r[,p.adjust_clusterProfiler:=NA_real_]
  cp <- suppressMessages(enricher(s,universe=u,TERM2GENE=as.data.frame(t2g),
    TERM2NAME=as.data.frame(unique(a[,.(ID,Description)])),pvalueCutoff=1,qvalueCutoff=1,minGSSize=10,maxGSSize=500,pAdjustMethod='BH'))
  if(!is.null(cp)) {
    q <- as.data.table(cp@result);ix <- match(q$ID,r$ID)
    stopifnot(!anyNA(ix),all(q$Count==r$Count[ix]),isTRUE(all.equal(q$pvalue,r$pvalue[ix],tolerance=1e-12,check.attributes=FALSE)))
    r[ix,p.adjust_clusterProfiler:=q$p.adjust]
  }
  r <- merge(r,unique(a[,.(ID,Description,Ontology)]),by='ID',all.x=TRUE)
  r[,significant:=FDR_fixed_family<.05];setorder(r,pvalue,ID);r
}

for(co in names(contrasts)) {
  v <- w[[contrasts[co]]];valid <- is.finite(v)&w$sites>=100
  if(grepl('^fd_',co))valid <- valid&v>0&v<1
  valid_sets[[co]] <- valid
  universe <- sort(unique(hits[window_id%in%w$window_id[valid],gene_id]))
  write_tsv(data.table(gene_id=universe),paste0(co,'_background_genes.tsv'))
  previous <- NULL
  for(p in ps) {
    cut <- as.numeric(quantile(v[valid],1-p/100,type=7));sel <- valid&v>=cut
    if(!is.null(previous))stopifnot(all(which(sel)%in%previous));previous <- which(sel)
    key <- paste(co,p,sep='_');selected_sets[[key]] <- sel
    gs <- sort(unique(hits[window_id%in%w$window_id[sel],gene_id]))
    write_tsv(w[sel,.(chr,start,end,window_id,statistic=get(contrasts[co]))],paste0(key,'_windows.tsv'))
    writeLines(gs,file.path(here,paste0(key,'_genes.txt')))
    s <- data.table(contrast=co,top_percent=p,cutoff=cut,eligible_windows=sum(valid),selected_windows=sum(sel),genes=length(gs))
    for(db in names(anno)) {
      er <- enrich(gs,universe,db);er[,`:=`(contrast=co,top_percent=p,database=db)]
      enrichments[[paste(key,db)]] <- er
      s[,(paste0(db,'_FDR05_terms')):=sum(er$significant)]
    }
    thresholds[[key]] <- s
    cat(key,':',sum(sel),'windows,',length(gs),'genes; GO significant',s$GO_FDR05_terms,'\n')
  }
}
thr <- rbindlist(thresholds);ers <- rbindlist(enrichments)
write_tsv(thr,'threshold_summary.tsv');write_tsv(ers,'enrichment_all.tsv')
write_tsv(ers[significant==TRUE],'enrichment_FDR05.tsv')
# Stratified genomic-block bootstrap; threshold labels fixed, no iid-window tests.
# Sample blocks within each chromosome, sharing the same bootstrap draws across
# traits and contrasts. Report percentile intervals for relative mean difference.
results <- list();B <- 2000L
for(block_bp in c(1000000L,5000000L)) {
  w[,block:=paste(chr,(start-1)%/%block_bp,sep=':')]
  btab <- unique(w[,.(block,chr)]);btab[,index:=.I]
  idx <- match(w$block,btab$block)
  weights <- matrix(0,nrow=B,ncol=nrow(btab))
  for(c in unique(btab$chr)) {
    bi <- btab[chr==c,index]
    weights[,bi] <- t(rmultinom(B,length(bi),rep(1/length(bi),length(bi))))
  }
  for(co in names(contrasts))for(p in ps)for(t in traits) {
    v <- w[[t]];eligible <- valid_sets[[co]]&is.finite(v)
    high <- selected_sets[[paste(co,p,sep='_')]]&eligible
    rest <- eligible&!high
    stopifnot(sum(high)>0,sum(rest)>0,mean(v[rest])>0)
    agg <- matrix(0,nrow=nrow(btab),ncol=4)
    for(j in 1:2) {
      mask <- if(j==1)high else rest
      temp <- data.table(id=idx[mask],value=v[mask])[,.(total=sum(value),n=.N),by=id]
      agg[temp$id,(j-1)*2+1] <- temp$total;agg[temp$id,(j-1)*2+2] <- temp$n
    }
    boot <- weights%*%agg
    effect <- 100*((boot[,1]/boot[,2])/(boot[,3]/boot[,4])-1)
    stopifnot(all(is.finite(effect)))
    ci <- quantile(effect,c(.025,.975),names=FALSE)
    results[[paste(block_bp,co,p,t)]] <- data.table(contrast=co,top_percent=p,trait=t,
      high_n=sum(high),background_n=sum(rest),mean_high=mean(v[high]),mean_background=mean(v[rest]),
      median_high=median(v[high]),median_background=median(v[rest]),
      relative_difference_percent=100*(mean(v[high])/mean(v[rest])-1),CI_low=ci[1],CI_high=ci[2],
      block_bp=block_bp,bootstrap_replicates=B,observed_blocks=nrow(btab))
  }
  cat('Completed',block_bp,'bp block bootstrap\n')
}
effects <- rbindlist(results);write_tsv(effects,'panel_c_effects_and_block_sensitivity.tsv')
write_tsv(w[,c('chr','start','end','sites','annotation_width','gene_coverage','repeat_coverage','NBS_density','NBS_count',unname(contrasts)),with=FALSE],
          'panel_c_input_windows.tsv')

writeLines(capture.output(sessionInfo()),file.path(here,'sessionInfo.txt'))
inputs <- c(src,annfile,nbsfile,file.path(root,'data/annotations/chromosome_lengths.tsv'),file.path(root,'data/annotations/gene_coordinates.tsv.gz'),
            file.path(prev,c('gene_GO.tsv.gz','gene_KEGG.tsv.gz')))
write_tsv(data.table(path=inputs,md5=unname(tools::md5sum(inputs))),'input_manifest.tsv')
# Rendering is a separate step after analysis completes.
cat('Analysis and panel rendering completed.\n')
