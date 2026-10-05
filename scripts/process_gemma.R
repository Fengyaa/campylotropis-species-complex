#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))
args <- commandArgs(trailingOnly=TRUE)
if(length(args)<3) stop('Usage: Rscript process_gemma.R TRAIT OUTDIR chr1.assoc.txt ... chr11.assoc.txt\nFor a premerged file, pass that one file. All 11 chromosomes must be represented.')
trait <- args[1]; out <- args[2]; files <- args[-c(1,2)]
if(anyDuplicated(normalizePath(files,mustWork=TRUE))) stop('Duplicated input files')
if(dir.exists(out)) stop('Output directory exists; use a new directory to avoid overwriting results')
required <- c('chr','rs','ps','n_miss','allele1','allele0','af','beta','se','p_wald')
parts <- lapply(files,function(f){
 d <- fread(f,na.strings=c('NA','NaN','nan','.'))
 if(!all(required %in% names(d))) stop(paste('Missing required columns:',f))
 # Reject multivariate outputs or malformed numbers instead of silently treating as univariate.
 for(k in c('ps','n_miss','af','beta','se','p_wald')) {
  old <- d[[k]]; x <- suppressWarnings(as.numeric(old))
  if(any(!is.na(old)&is.na(x))) stop(paste('Nonnumeric column',k,'in',f))
  set(d,j=k,value=x)
 }
 if('pve' %in% names(d)) d[,pve:=NULL]
 d[,source_file:=basename(f)];d
})
d <- rbindlist(parts,use.names=TRUE,fill=TRUE);rm(parts);gc(verbose=FALSE)
if(nrow(d)==0) stop('No tests')
if(anyNA(d$rs)||anyNA(d$ps)) stop('Missing SNP ID or position')
# IDs may be scaffold:position or scaffold:position:REF:ALT.
parsed <- regexec('^(HiC_scaffold([0-9]+)):([0-9]+)(:.*)?$',d$rs)
match <- regmatches(d$rs,parsed)
if(any(lengths(match)==0)) stop('Unexpected rs format: expected HiC_scaffoldN:position')
d[,scaffold:=vapply(match,`[`,character(1),2)]
d[,chromosome:=as.integer(vapply(match,`[`,character(1),3))]
idpos <- as.numeric(vapply(match,`[`,character(1),4));rm(parsed,match)
if(any(idpos!=d$ps)) stop('rs and ps disagree')
if(!setequal(unique(d$chromosome),1:11)) stop('Not all 11 chromosomes are represented; do not apply per-chromosome thresholds')
if(anyDuplicated(d$rs)) stop('Repeated rs IDs; check duplicated files/records before processing')
if(anyDuplicated(d[,.(scaffold,ps)])) stop('Repeated positions; confirm biallelic test design before processing')
if(any(d$p_wald<0|d$p_wald>1,na.rm=TRUE)) stop('P outside [0,1]')
M <- nrow(d); threshold <- .05/M; suggestive <- 1/M
# Count all reported tests, including failures, in Bonferroni denominator.
d[,valid_test:=is.finite(p_wald)&is.finite(beta)&is.finite(se)&se>0]
d[,P_BONF:=NA_real_];d[valid_test==TRUE,P_BONF:=pmin(1,p_wald*M)]
d[,Q_BH:=NA_real_];d[valid_test==TRUE,Q_BH:=p.adjust(p_wald,method='BH',n=M)]
d[,significant_bonferroni:=valid_test & !is.na(p_wald) & p_wald<=threshold]
d[,suggestive_1_over_M:=valid_test & !is.na(p_wald) & p_wald<=suggestive]
d[,zero_p_flag:=!is.na(p_wald)&p_wald==0]
setnames(d,c('allele1','allele0','af'),c('effect_allele','other_allele','effect_allele_frequency'))
setorder(d,chromosome,ps)
dir.create(out,recursive=TRUE)
fwrite(d,file.path(out,paste0(trait,'.all_tests.tsv.gz')),sep='\t',na='NA')
fwrite(d[significant_bonferroni==TRUE],file.path(out,paste0(trait,'.bonferroni.tsv')),sep='\t',na='NA')
fwrite(d[suggestive_1_over_M==TRUE],file.path(out,paste0(trait,'.suggestive.tsv')),sep='\t',na='NA')
fwrite(d[valid_test==FALSE],file.path(out,paste0(trait,'.failed_tests.tsv')),sep='\t',na='NA')
top <- d[valid_test==TRUE][order(p_wald,chromosome,ps)][seq_len(min(100,sum(d$valid_test)))]
fwrite(top,file.path(out,paste0(trait,'.top100.tsv')),sep='\t',na='NA')
summary <- data.table(trait=trait,input_files=length(files),reported_tests=M,valid_tests=sum(d$valid_test),failed_tests=sum(!d$valid_test),threshold_bonferroni=threshold,threshold_suggestive=suggestive,significant_bonferroni=sum(d$significant_bonferroni),suggestive=sum(d$suggestive_1_over_M),zero_P=sum(d$zero_p_flag))
fwrite(summary,file.path(out,'summary.tsv'),sep='\t');print(summary)
fwrite(d[,.(reported_tests=.N,valid_tests=sum(valid_test),significant=sum(significant_bonferroni)),by=chromosome],file.path(out,'chromosome_summary.tsv'),sep='\t')
writeLines(c('Single-trait genome-wide correction; not a joint correction across eight traits.','1/M is suggestive, not FWER=0.05.','effect_allele/other_allele are GEMMA allele1/allele0, not genomic REF/ALT.','No SNP PVE calculated. Existing added pve column discarded.','No observed sample count inferred from n_miss: confirm analyzed sample count in GEMMA logs.','Original P values kept, including zero underflow flags.','No causal interpretation or independent locus count inferred from significant SNP count.',paste('Inputs:',paste(normalizePath(files),collapse='; '))),file.path(out,'notes.txt'))
