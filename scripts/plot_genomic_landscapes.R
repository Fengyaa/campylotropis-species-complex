# Paths are relative to this archive, independent of the caller's working directory.
.script <- sub('^--file=', '', commandArgs(FALSE)[grepl('^--file=', commandArgs(FALSE))])
root <- normalizePath(file.path(dirname(.script), '..'))
setwd(root)
here <- 'results/windows';dir.create(here,recursive=TRUE,showWarnings=FALSE)
suppressPackageStartupMessages(library(data.table))
src <- 'data/windows/genomic_windows.tsv.gz'
gp <- 'data/annotations/gene_coordinates.tsv.gz'
np <- 'data/annotations/NBS_gene_ids.txt'
d <- fread(src,na.strings=c('NA','nan','NaN',''))
g <- fread(gp,header=FALSE,col.names=c('chr','gene','start','end'))
nbs <- unique(sub('[.]m[0-9]+$','',trimws(readLines(np))))
stopifnot(all(nbs%in%g$gene),!anyDuplicated(g$gene))
d[,chromosome:=as.integer(sub('HiC_scaffold','',chr))]
setorder(d,chromosome,start)
ends <- d[,.(endpoint=max(end)),by=.(chr,chromosome)]
ends <- merge(ends,fread('data/annotations/chromosome_lengths.tsv'),by='chr')
ends[,endpoint:=reference_length]
setorder(ends,chromosome)
# Fixed 100-kb annotation bins, independent of availability of FST/fd estimates.
bins <- rbindlist(lapply(1:nrow(ends),function(i) {
  ss <- seq(1,ends$endpoint[i],by=100000)
  data.table(chr=ends$chr[i],start=ss,end=pmin(ss+99999,ends$endpoint[i]))
}))
bins[,bin_id:=.I]
# Union gene intervals before measuring coverage: overlapping genes counted once.
setorder(g,chr,start,end)
g[,group:=cumsum(start>shift(cummax(as.numeric(end)),fill=-Inf)+1),by=chr]
unions <- g[,.(start=min(start),end=max(end)),by=.(chr,group)][,group:=NULL]
setkey(unions,chr,start,end)
h <- foverlaps(bins,unions,type='any',nomatch=0L)
cov <- h[,.(gene_bp=sum(pmin(end,i.end)-pmax(start,i.start)+1)),by=bin_id]
bins <- merge(bins,cov,by='bin_id',all.x=TRUE)
bins[is.na(gene_bp),gene_bp:=0]
ng <- copy(g[gene%in%nbs,.(chr,start,end,gene)])
setkey(ng,chr,start,end)
nh <- foverlaps(bins[,.(chr,start,end,bin_id)],ng,by.x=c('chr','start','end'),by.y=c('chr','start','end'),type='any',nomatch=0L)
nc <- nh[,.(nbs_count=uniqueN(gene)),by=bin_id]
bins <- merge(bins,nc,by='bin_id',all.x=TRUE)
bins[is.na(nbs_count),nbs_count:=0]
bins[,gene_fraction:=gene_bp/(end-start+1)]
stopifnot(all(bins$gene_fraction>=0 & bins$gene_fraction<=1),all(bins$nbs_count>=0))
fwrite(bins,file.path(here,'annotation_gene_NBS_matched_100kb.tsv'),sep='\t')
fwrite(g[gene%in%nbs,.(chr,gene,start,end)],file.path(here,'NBS_gene_coordinates.tsv'),sep='\t')
pal_gene <- colorRampPalette(c('#F4F5E9','#9CCDB5','#267C83','#164451'))(101)
pal_nbs <- colorRampPalette(c('#F3F1F6','#B1A0CC','#63358C','#31114E'))(101)
nmax <- max(bins$nbs_count)
cols <- c('#315F78','#BD673B')
limits_for <- function(kind) list(c(0,if(kind=='species') .3 else .5),c(0,1))
# Clamp negative estimates only for display; preserve input and missing values.
plot_value <- function(v) pmax(v, 0)
all_fields <- c('Fst_Cpo_Cyu','species_trio_fd_recomputed','Mean_fst_sym','Mean_f_d_sym_2026')
audit <- rbindlist(lapply(seq_along(all_fields),function(i) {
  v <- d[[all_fields[i]]]; z <- plot_value(v)
  lim <- limits_for(if(i<=2) 'species' else 'sympatric')[[if(i %% 2 == 1) 1 else 2]]
  stopifnot(all(z[is.finite(z)] >= lim[1]),all(z[is.finite(z)] <= lim[2]),
            identical(is.na(v),is.na(z)))
  data.table(field=all_fields[i],finite_windows=sum(is.finite(v)),
             negative_windows_set_to_zero=sum(v<0,na.rm=TRUE),
             axis_min=lim[1],axis_max=lim[2],missing_windows=sum(is.na(v)),display_min=min(z,na.rm=TRUE),display_max=max(z,na.rm=TRUE))
}))
fwrite(audit,file.path(here,'negative_value_display_audit_180mm.tsv'),sep='\t')
# Single drawing canvas keeps 44 tracks compact and axes explicitly aligned.
draw <- function(kind) {
  limits <- limits_for(kind)
  fields <- if(kind=='species')c('Fst_Cpo_Cyu','species_trio_fd_recomputed') else c('Mean_fst_sym','Mean_f_d_sym_2026')
  par(mar=c(0,0,0,0),family='Helvetica',xpd=NA)
  plot.new();plot.window(xlim=c(-13,123),ylim=c(-.2,13.3),xaxs='i',yaxs='i')
  text(-12.2,13.03,if(kind=='species')'Species-level differentiation and introgression' else 'Sympatric differentiation and introgression',adj=0,cex=.9,font=2)
  text(-12.2,12.73,if(kind=='species')'FST: Cpo-Cyu; fd: Cma-Cyu-Cpo trio | 10-kb summaries' else 'Mean sympatric FST and corrected mean fd | 10-kb summaries',adj=0,cex=.65,col='#53616B')
  ticks <- seq(0,120,20)
  for(i in 1:11) {
    c <- paste0('HiC_scaffold',i);a <- d[chr==c];b <- bins[chr==c]
    e <- ends[chr==c,endpoint]/1e6;base <- 12.15-(i-1)*1.035
    text(-5.8,base-.42,paste0('Chr',i),adj=1,font=2,cex=.7)
    # Annotation strips identified by the shared legend above each pair of continuous tracks.
    for(k in 1:2) {
      yl <- base-(k-1)*.105
      val <- list(b$gene_fraction,b$nbs_count/nmax)[[k]]
      pal <- list(pal_gene,pal_nbs)[[k]]
      fill <- pal[pmin(101,floor(val*100)+1)]
      fill[is.na(fill)] <- '#AEB4BA'
      rect((b$start-1)/1e6,yl-.075,b$end/1e6,yl,col=fill,border=NA)
    }
    for(j in 1:2) {
      top <- base-.225-(j-1)*.33;bottom <- top-.285
      lim <- limits[[j]]
      yy <- function(v)bottom+(v-lim[1])/diff(lim)*.285
      rect(0,bottom,e,top,col=if(j==1)'#F2F6F8' else '#FBF4EE',border=NA)
      segments(ticks[ticks<=e],bottom,ticks[ticks<=e],top,col='#E3E7E9',lwd=.35)
      segments(0,yy(0),e,yy(0),col='#A0A9AE',lwd=.35)
      v <- plot_value(a[[fields[j]]]); xx <- (a$start+a$end)/2e6
      gap <- c(FALSE,a$start[-1]>head(a$end,-1)+1)
      ii <- rep(seq_along(xx),1+as.integer(gap));xp <- xx[ii];yp <- yy(v[ii])
      first <- !duplicated(ii)&gap[ii];xp[first]<-NA_real_;yp[first]<-NA_real_
      lines(xp,yp,col=cols[j],lwd=.25)
      points(xp,yp,pch=16,cex=.035,col=cols[j])
      text(-.8,(top+bottom)/2,if(j==1)expression(F[ST]) else expression(f[d]),adj=1,cex=.6,col=cols[j])
    }
  }
  segments(0,12.32,120,12.32,col='#AAB2B8',lwd=.5)
  segments(ticks,12.32,ticks,12.38,col='#AAB2B8',lwd=.5)
  text(ticks,12.48,ticks,cex=.58,col='#46515A')
  segments(0,.83,120,.83,col='#AAB2B8',lwd=.5)
  segments(ticks,.83,ticks,.77,col='#AAB2B8',lwd=.5)
  text(ticks,.62,ticks,cex=.58,col='#46515A')
  text(60,.34,'Chromosome position (Mb)',cex=.7)
  # Legends occupy the whitespace to the right of shorter chromosomes.
  bar <- function(x,y,pal,title,labels) {
    text(x,y+.26,title,adj=0,cex=.58,col='#34424B')
    xs <- seq(x,x+25,length.out=102)
    rect(xs[-102],y,xs[-1],y+.1,col=pal,border=NA)
    text(c(x,x+12.5,x+25),y-.15,labels,cex=.5,col='#52616B')
  }
  bar(82,5.25,pal_gene,'Gene-body coverage',c('0','0.5','1'))
  bar(82,4.48,pal_nbs,'NBS genes / bin',c('0',format(nmax/2,trim=TRUE),as.character(nmax)))
  text(82,3.84,'Annotation bins: 100 kb',adj=0,cex=.52,col='#52616B')
  text(82,3.57,'Statistics: 10 kb',adj=0,cex=.52,col='#52616B')
  text(82,3.30,paste0('FST: 0 to ',limits[[1]][2],'; fd: 0 to 1'),adj=0,cex=.5,col='#52616B')
  text(82,3.03,'Negative estimates plotted as 0',adj=0,cex=.48,col='#52616B')
  text(-12.2,-.04,'Gene = union gene-body coverage; NBS = existing BLAST candidate list. Missing statistics remain gaps.',adj=0,cex=.46,col='#52616B')
}
for(kind in c('species','sympatric')) {
  stem <- file.path(here,paste0('Fig3_',kind,'_matched_gene_NBS_180mm'))
  pdf(paste0(stem,'.pdf'),width=180/25.4,height=140/25.4,useDingbats=FALSE);draw(kind);dev.off()
  png(paste0(stem,'.png'),width=180,height=140,units='mm',res=260,type='cairo');draw(kind);dev.off()
}
fwrite(data.table(path=c(src,gp,np),md5=unname(tools::md5sum(c(src,gp,np)))),
       file.path(here,'gene_NBS_matched_input_manifest.tsv'),sep='\t')
cat('Saved both matched Gene/NBS figures; NBS genes:',length(nbs),'\n')
