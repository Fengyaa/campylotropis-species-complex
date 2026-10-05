# Paths are relative to this archive, independent of the caller's working directory.
.script <- sub('^--file=', '', commandArgs(FALSE)[grepl('^--file=', commandArgs(FALSE))])
root <- normalizePath(file.path(dirname(.script), '..'))
setwd(root)
suppressPackageStartupMessages({library(data.table);library(ggplot2);library(grid)})
here <- 'results/windows'
e <- fread(file.path(here,'panel_c_effects_and_block_sensitivity.tsv'))
a <- fread(file.path(here,'enrichment_all.tsv'))
cols <- c(FST='#315F78',fd='#BD673B')
style <- theme_classic(base_size=9,base_family='Helvetica')+
 theme(plot.title=element_text(size=11,face='bold'),plot.subtitle=element_text(size=8,color='#596873'),
 plot.title.position='plot',legend.position='bottom',strip.background=element_blank(),strip.text=element_text(size=8,face='bold'),
 axis.text=element_text(color='#42515C'),plot.margin=margin(5,7,5,5),panel.grid.major.y=element_line(color='#E8ECEF',linewidth=.25))
e[,`:=`(mean_ratio=1+relative_difference_percent/100,ratio_CI_low=1+CI_low/100,ratio_CI_high=1+CI_high/100)]
stopifnot(isTRUE(all.equal(e$mean_ratio,e$mean_high/e$mean_background,tolerance=1e-10)))
fwrite(e,file.path(here,'panel_c_effects_as_ratios.tsv'),sep='\t')
e <- e[block_bp==1000000]
e[,context:=factor(ifelse(grepl('sympatric',contrast),'Sympatric','Species'),levels=c('Species','Sympatric'))]
e[,metric:=factor(ifelse(grepl('^FST',contrast),'FST','fd'),levels=c('FST','fd'))]
e[,threshold:=factor(paste0(top_percent,'%'),levels=c('10%','5%','2.5%','1%'))]
e[,feature:=factor(trait,levels=c('gene_coverage','repeat_coverage','NBS_density'),labels=c('Gene-body coverage','Repeat coverage','NBS genes per Mb'))]
e[,panel:=factor(paste(context,feature,sep=' | '),levels=as.vector(t(outer(c('Species','Sympatric'),c('Gene-body coverage','Repeat coverage','NBS genes per Mb'),paste,sep=' | '))))]
dodge <- position_dodge(.25)
pc <- ggplot(e,aes(threshold,mean_ratio,color=metric,group=metric))+
 geom_hline(yintercept=1,linetype='dashed',color='#9EAAB1',linewidth=.35)+
 geom_line(position=dodge,linewidth=.4)+geom_errorbar(aes(ymin=ratio_CI_low,ymax=ratio_CI_high),position=dodge,width=.13,linewidth=.35)+geom_point(position=dodge,size=1.65)+
 facet_wrap(~panel,ncol=3,scales='free_y')+scale_color_manual(values=cols,breaks=names(cols),labels=c('High FST','High fd'),name=NULL)+
 labs(title='c  Genomic features of candidate regions',subtitle='Dashed line = background (1-fold); 95% block-bootstrap intervals; panel y scales differ',
 x='Upper-tail candidate threshold',y='Candidate / background mean ratio')+style
# Show the union of terms significant in either context for each statistic.
a <- a[database=='GO']
a[,metric:=ifelse(grepl('^FST',contrast),'FST','fd')]
a[,context:=factor(ifelse(grepl('sympatric',contrast),'Sympatric','Species'),levels=c('Species','Sympatric'))]
t <- a[significant==TRUE,.(best=min(FDR_fixed_family)),by=.(metric,ID,Description)]
setorder(t,metric,best)
z <- merge(a,t[,.(metric,ID)],by=c('metric','ID'))
stopifnot(all(z[,.N,by=.(metric,ID)]$N==8))
wrap <- function(s)paste(strwrap(s,width=31),collapse='\n')
lab <- vapply(t$Description,wrap,character(1))
# Make each metric/term key unique, including shared GO IDs if present.
t[,key:=paste(metric,ID,sep=':')];z[,key:=paste(metric,ID,sep=':')]
z[,term:=factor(key,levels=rev(t$key),labels=rev(lab))]
z[,metric_label:=factor(metric,levels=c('FST','fd'),labels=c('High FST','High fd'))]
z[,threshold:=factor(paste0(top_percent,'%'),levels=c('10%','5%','2.5%','1%'))]
z[,log2fold:=log2(pmax(FoldEnrichment,.25))]
z[,sig:=factor(ifelse(significant,'FDR < 0.05','Not significant'),levels=c('Not significant','FDR < 0.05'))]
fwrite(z,file.path(here,'panel_d_displayed_terms.tsv'),sep='\t')
mx <- max(z$log2fold)
pd <- ggplot(z,aes(threshold,term))+
 geom_point(aes(size=Count,fill=log2fold,shape=sig),color='#3E4952',stroke=.35)+
 facet_grid(rows=vars(metric_label),cols=vars(context),scales='free_y',space='free_y')+
 scale_shape_manual(values=c('Not significant'=21,'FDR < 0.05'=23),name=NULL)+
 scale_fill_gradientn(colors=c('#D5E2EE','#FAFAFA','#F1BC84','#A54316'),values=scales::rescale(c(-2,0,2,mx),from=c(-2,mx)),
 limits=c(-2,mx),breaks=c(0,1,2,3,4)[c(0,1,2,3,4)<=mx],labels=c('1','2','4','8','16')[c(0,1,2,3,4)<=mx],name='Enrichment fold')+
 scale_size_area(max_size=5,breaks=c(5,50,200),name='Hit genes')+
 labs(title='d  Functional enrichment across thresholds',subtitle='All GO terms significant at any threshold in either context; fixed-family BH correction',
 x='Upper-tail candidate threshold',y=NULL)+style+
 theme(axis.text.y=element_text(size=7),strip.text.y=element_text(angle=0,size=8),axis.line.y=element_blank(),axis.ticks.y=element_blank(),
 panel.grid.major.y=element_line(color='#F0F1F2',linewidth=.2),legend.text=element_text(size=7),legend.title=element_text(size=7),panel.spacing=unit(4,'mm'))+
 guides(fill=guide_colorbar(order=1,barwidth=unit(25,'mm'),barheight=unit(2,'mm'),title.position='top'),
 size=guide_legend(order=2,title.position='top'),shape=guide_legend(order=3,ncol=1))
savefig <- function(p,stem,h){ggsave(file.path(here,paste0(stem,'.pdf')),p,width=190,height=h,units='mm',useDingbats=FALSE);ggsave(file.path(here,paste0(stem,'.png')),p,width=190,height=h,units='mm',dpi=260,bg='white')}

# Publication panels b and c, with species panels preceding sympatric panels.
e[,panel:=factor(paste(context,feature,sep=' | '),levels=as.vector(t(outer(c('Species','Sympatric'),c('Gene-body coverage','Repeat coverage','NBS genes per Mb'),paste,sep=' | '))))]
pb <- ggplot(e,aes(threshold,mean_ratio,color=metric,group=metric))+
 geom_hline(yintercept=1,linetype='dashed',color='#9EAAB1',linewidth=.3)+
 geom_line(position=dodge,linewidth=.35)+geom_errorbar(aes(ymin=ratio_CI_low,ymax=ratio_CI_high),position=dodge,width=.1,linewidth=.3)+geom_point(position=dodge,size=1.2)+
 facet_wrap(~panel,nrow=1,scales='free_y',labeller=as_labeller(c('Species | Gene-body coverage'='Species\nGene coverage','Species | Repeat coverage'='Species\nRepeat coverage','Species | NBS genes per Mb'='Species\nNBS density','Sympatric | Gene-body coverage'='Sympatric\nGene coverage','Sympatric | Repeat coverage'='Sympatric\nRepeat coverage','Sympatric | NBS genes per Mb'='Sympatric\nNBS density')))+scale_color_manual(values=cols,labels=c('High FST','High fd'),name=NULL)+labs(title='b',x=NULL,y='Candidate / background')+style+
 theme(plot.title=element_text(size=9),strip.text=element_text(size=5.6),axis.text.x=element_text(size=5.5,angle=45,hjust=1),axis.text.y=element_text(size=5.5),axis.title.y=element_text(size=6.5),panel.spacing=unit(1,'mm'),legend.text=element_text(size=7),plot.margin=margin(2,2,1,2))
pgo <- pd+labs(title='c',subtitle=NULL,x=NULL)+theme(plot.title=element_text(size=9),axis.text.y=element_text(size=6.5),strip.text=element_text(size=7),legend.text=element_text(size=6.5))
ggsave(file.path(here,'Figure_3b.pdf'),pb,width=180,height=49,units='mm',useDingbats=FALSE)
ggsave(file.path(here,'Figure_3c.pdf'),pgo,width=180,height=max(70,40+nrow(t)*6),units='mm',useDingbats=FALSE)
