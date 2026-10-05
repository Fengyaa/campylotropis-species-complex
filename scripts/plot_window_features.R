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
savefig(pc,'Fig3c_species_sympatric_features',105)
savefig(pd,'Fig3d_species_sympatric_GO',145)
draw <- function(){grid.newpage();print(pc,vp=viewport(x=.5,y=.79,width=1,height=.42));print(pd,vp=viewport(x=.5,y=.29,width=1,height=.58))}
pdf(file.path(here,'Fig3cd_revised.pdf'),width=190/25.4,height=250/25.4,useDingbats=FALSE);draw();dev.off()
png(file.path(here,'Fig3cd_revised.png'),width=190,height=250,units='mm',res=260,type='cairo');draw();dev.off()
# Manuscript-width compact layout.  Keep the full-size panels above for inspection.
pc_compact <- ggplot(e,aes(threshold,mean_ratio,color=metric,group=metric))+
 geom_hline(yintercept=1,linetype='dashed',color='#9EAAB1',linewidth=.3)+
 geom_line(position=dodge,linewidth=.35)+
 geom_errorbar(aes(ymin=ratio_CI_low,ymax=ratio_CI_high),position=dodge,width=.11,linewidth=.3)+
 geom_point(position=dodge,size=1.2)+
 facet_grid(rows=vars(feature),cols=vars(context),scales='free_y')+
 scale_color_manual(values=cols,breaks=names(cols),labels=c('High FST','High fd'),name=NULL)+
 labs(title='c  Candidate-region features',x=NULL,y='Candidate / background')+style+
 theme(legend.position='none',plot.title=element_text(size=8.5),strip.text=element_text(size=6.5),
 strip.text.y=element_text(angle=0),axis.text=element_text(size=5.8),axis.title.y=element_text(size=6.5),
 panel.spacing=unit(1.2,'mm'),plot.margin=margin(2,2,0,2))
pd_compact <- pd + labs(title='d  GO enrichment',subtitle=NULL,x=NULL)+
 scale_size_area(max_size=3.2,breaks=c(5,50,200),name='Hit genes')+
 theme(legend.position='none',plot.title=element_text(size=8.5),strip.text=element_text(size=6.5),
 axis.text.x=element_text(size=6),axis.text.y=element_text(size=5.7,lineheight=.85),
 panel.spacing=unit(1.2,'mm'),plot.margin=margin(0,2,2,2))
draw_compact <- function(){
 grid.newpage()
 print(pc_compact,vp=viewport(x=.5,y=1-40/98/2,width=1,height=40/98))
 print(pd_compact,vp=viewport(x=.5,y=55/98/2,width=1,height=55/98))
}
pdf(file.path(here,'Fig3cd_compact_98mm.pdf'),width=190/25.4,height=98/25.4,useDingbats=FALSE);draw_compact();dev.off()
png(file.path(here,'Fig3cd_compact_98mm.png'),width=190,height=98,units='mm',res=300,type='cairo');draw_compact();dev.off()
# Alternative requested layout: all six feature panels in one horizontal row.
e[,panel_one_row:=factor(paste(context,feature,sep=' | '),
 levels=as.vector(outer(c('Species','Sympatric'),c('Gene-body coverage','Repeat coverage','NBS genes per Mb'),paste,sep=' | ')),
 labels=c('Species · Gene','Sympatric · Gene','Species · Repeat','Sympatric · Repeat','Species · NBS','Sympatric · NBS'))]
pc_one_row <- ggplot(e,aes(threshold,mean_ratio,color=metric,group=metric))+
 geom_hline(yintercept=1,linetype='dashed',color='#9EAAB1',linewidth=.3)+
 geom_line(position=dodge,linewidth=.35)+
 geom_errorbar(aes(ymin=ratio_CI_low,ymax=ratio_CI_high),position=dodge,width=.1,linewidth=.3)+
 geom_point(position=dodge,size=1.2)+
 facet_wrap(~panel_one_row,nrow=1,scales='free_y')+
 scale_color_manual(values=cols,breaks=names(cols),labels=c('High FST','High fd'),name=NULL)+
 labs(title='c',x=NULL,y='Candidate / background')+style+
 theme(legend.position='none',plot.title=element_text(size=8.5),strip.text=element_text(size=5.7),
 axis.text.x=element_text(size=5.3,angle=45,hjust=1),axis.text.y=element_text(size=5.2),
 axis.title.y=element_text(size=6.3),panel.spacing=unit(.8,'mm'),plot.margin=margin(2,2,0,2))
draw_one_row <- function(){
 grid.newpage()
 print(pc_one_row,vp=viewport(x=.5,y=1-29/98/2,width=1,height=29/98))
 print(pd_compact,vp=viewport(x=.5,y=66/98/2,width=1,height=66/98))
}
pdf(file.path(here,'Fig3cd_6panels_one_row_98mm.pdf'),width=190/25.4,height=98/25.4,useDingbats=FALSE);draw_one_row();dev.off()
png(file.path(here,'Fig3cd_6panels_one_row_98mm.png'),width=190,height=98,units='mm',res=300,type='cairo');draw_one_row();dev.off()
pc_one_row_standalone <- pc_one_row+
 labs(title='c  Candidate-region features')+
 theme(legend.position='bottom',legend.direction='horizontal',legend.text=element_text(size=6.5),
 legend.key.width=unit(6,'mm'),legend.margin=margin(0,0,0,0),
 plot.margin=margin(2,2,1,2))
ggsave(file.path(here,'Fig3c_6panels_one_row_180mm.pdf'),pc_one_row_standalone,
 width=180,height=46,units='mm',useDingbats=FALSE)
ggsave(file.path(here,'Fig3c_6panels_one_row_180mm.png'),pc_one_row_standalone,
 width=180,height=46,units='mm',dpi=300,bg='white')
ggsave(file.path(here,'Fig3c_6panels_one_row_170mm.pdf'),pc_one_row_standalone,
 width=170,height=46,units='mm',useDingbats=FALSE)
ggsave(file.path(here,'Fig3c_6panels_one_row_170mm.png'),pc_one_row_standalone,
 width=170,height=46,units='mm',dpi=300,bg='white')
cat('Displayed',nrow(t),'GO terms across both contexts.\n')
