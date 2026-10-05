suppressPackageStartupMessages(library(data.table))
.script <- sub('^--file=', '', commandArgs(FALSE)[grepl('^--file=', commandArgs(FALSE))])
root <- normalizePath(file.path(dirname(.script), '..'));setwd(root)

out<-'results/gea'
r<-fread(file.path(out,'comparison_statistics.tsv'))[candidate_set=='Z3.5_r0.5']
ci<-fread(file.path(out,'block_bootstrap_intervals.tsv'))[candidate_set=='Z3.5_r0.5'&block_bp==1e6]
traits<-c(setdiff(r$trait,c('Rho','species_trio_fd_recomputed')),'species_trio_fd_recomputed')
labels<-c('Diversity: CMA','Diversity: CPO','Diversity: CYU','FST: CMA-CPO','FST: CMA-CYU','FST: CPO-CYU','DXY: CMA-CPO','DXY: CMA-CYU','DXY: CPO-CYU','XP-CLR: CMA / CPOCYU','XP-CLR: CPO / CMACYU','XP-CLR: CYU / CMACPO','XP-CLR: CPO / CYU','XP-CLR: CYU / CPO','Diversity: CPO allopatric','Diversity: CPO sympatric','Diversity: CYU sympatric','FST: sympatric mean','DXY: sympatric mean','fd: sympatric mean','fd: CPO-CYU')
a<-merge(r,ci,by=c('candidate_set','trait'));a<-a[match(traits,trait)];stopifnot(nrow(a)==21,tail(a$trait,1)=='species_trio_fd_recomputed',!any(a$trait=='Rho'))
fwrite(a,file.path(out,'Z3.5_r0.5_90mm_plot_data.tsv'),sep='\t')
fun<-function(){
 par(family='Helvetica',ps=7,mai=c(.64,1.46,.47,.10),mgp=c(1.7,.4,0),tcl=-.18,xaxs='i')
 yy<-21:1;xr<-range(c(a$standardized_ci_low,a$standardized_ci_high,0))+c(-.05,.05)
 plot(0,0,type='n',xlim=xr,ylim=c(.5,21.5),axes=FALSE,xlab='',ylab='')
 for(y in seq(2,20,2))rect(xr[1],y-.5,xr[2],y+.5,col='#F3F5F7',border=NA)
 abline(v=0,lty=2,lwd=.65,col='#8D969E')
 cols<-ifelse(a$q_BH_within_set<.05,'#B34D32','#647989')
 segments(a$standardized_ci_low,yy,a$standardized_ci_high,yy,col=cols,lwd=.8)
 points(a$standardized_difference,yy,col=cols,pch=19,cex=.65)
 axis(2,at=yy,labels=labels,las=2,tick=FALSE,cex.axis=.83)
 axis(1,at=pretty(xr,4),cex.axis=.88);box(col='#AAB1B7',lwd=.6)
 mtext('Standardized mean difference',side=1,line=1.45,cex=.88)
 mtext('Z > 3.5; |r| > 0.5',side=3,line=1.5,font=2,cex=1.1)
 mtext(paste0('pRDA-BayPass matched windows (n = ',format(fread(file.path(out,'window_qc.tsv'))[candidate_set=='Z3.5_r0.5',present_in_parameter_table],big.mark=','),')'),side=3,line=.35,cex=.81)
 par(xpd=NA)
 legend('bottomleft',inset=c(-.55,-.13),legend=c('BH-adjusted P < 0.05','Not significant'),col=c('#B34D32','#647989'),pch=19,pt.cex=.7,cex=.8,bty='n',ncol=2,x.intersp=.5)
}
pdf(file.path(out,'Z3.5_r0.5_parameters_90mm.pdf'),width=90/25.4,height=160/25.4,pointsize=7,useDingbats=FALSE);fun();dev.off()
png(file.path(out,'Z3.5_r0.5_parameters_90mm.png'),width=90,height=160,units='mm',res=600,pointsize=7,type='cairo');fun();dev.off()
writeLines(c('Width: 90 mm; height: 160 mm. Six-climate-variable pRDA Z>3.5, absolute Pearson r>0.5, shared with BayPass.', 'Recombination omitted; species-trio fd last. Values and significance unchanged: BH remains across all original 22 tests, not recalculated for 21 plotted rows.', 'Points: candidate-minus-background mean difference divided by genome-wide SD. Bars: chromosome-stratified 1-Mb block-bootstrap 95% intervals.', 'Red: 9,999-shift spatial permutation P adjusted by BH within candidate set <0.05. Spatial P and bootstrap intervals answer different questions.'),file.path(out,'Z3.5_r0.5_parameters_90mm_caption.txt'))
