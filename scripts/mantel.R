# Paths are relative to this archive, independent of the caller's working directory.
.script <- sub('^--file=', '', commandArgs(FALSE)[grepl('^--file=', commandArgs(FALSE))])
root <- normalizePath(file.path(dirname(.script), '..'))
setwd(root)
#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(vegan))
base <- 'data/spatial'
out <- 'results/mantel';dir.create(out,recursive=TRUE,showWarnings=FALSE)
seed <- 20261001L
nperm <- 9999L
fst_path <- file.path(base,'fst_matrix_WC84_diploid.csv')
fst <- as.matrix(read.csv(fst_path,row.names=1,check.names=FALSE))
ids <- rownames(fst)
stopifnot(identical(ids,colnames(fst)),isSymmetric(fst),all(is.finite(fst)),all(fst<1))
linear <- fst/(1-fst)
diag(linear) <- 0
write.csv(linear,file.path(out,'fst_matrix_linearized.csv'),quote=FALSE)
# The old distance files have numerical row labels. Recover population labels
# from their original 33-population ordering, independently checked below.
coords <- read.delim(file.path(base,'pop_coord.txt'),check.names=FALSE)
canonical <- function(x) gsub('-','',x,fixed=TRUE)
stopifnot(identical(coords$order,seq_len(nrow(coords))),
          identical(canonical(coords$pop),canonical(ids)))
clim <- read.csv(file.path(base,'coords_climate_all.csv'),check.names=FALSE)
stopifnot(identical(canonical(clim$code),canonical(ids)),
          max(abs(clim$Longitude-coords$Longitude))<1e-6,
          max(abs(clim$Latitude-coords$Latitude))<1e-6)
files <- c('Geographic distance'='allpop_distance_matrix.txt',
           'Environmental distance'='allpop_allenv_PCA_matrix.txt',
           'Precipitation'='allpop_precip_PCA_matrix.txt',
           'Temperature'='allpop_temp_PCA_matrix.txt',
           'Solar radiation'='allpop_rsrad_PCA_matrix.txt',
           'Wind speed'='allpop_rwind_PCA_matrix.txt',
           'Water vapor pressure'='allpop_rvapr_PCA_matrix.txt')
distances <- lapply(files,function(f) {
  m <- as.matrix(read.table(file.path(base,f),header=TRUE,check.names=FALSE))
  stopifnot(identical(dim(m),dim(fst)),all(is.finite(m)),
            max(abs(m-t(m)))<1e-8,max(abs(diag(m)))<1e-8)
  dimnames(m) <- list(ids,ids)
  m
})
# Preserve the active (not commented out) indices in mantel_test.R.
groups <- list(Total=seq_along(ids),Cma=c(4,5,23,24,25,26,27,28,33),
               Cyu=c(6,10,15,17,19,21,22,31,30),
               Cpo=c(1,3,7,8,9,12,14,16,20,29,32))
membership <- data.frame(order=seq_along(ids),population=ids,
                         Longitude=coords$Longitude,Latitude=coords$Latitude,
                         species_test='Total only')
for(g in names(groups)[-1]) membership$species_test[groups[[g]]] <- g
write.csv(membership,file.path(out,'population_groups.csv'),row.names=FALSE,quote=FALSE)
all_res <- list(); perm_sets <- list(); k <- 0L
for(gi in seq_along(groups)) {
 g <- names(groups)[gi]; ix <- groups[[gi]]
 set.seed(seed+gi-1L)
 perm <- permute::shuffleSet(length(ix),nset=nperm)
 perm_sets[[g]] <- perm
 y <- as.dist(linear[ix,ix]); geo <- as.dist(distances[[1]][ix,ix])
 for(j in seq_along(distances)) {
  x <- as.dist(distances[[j]][ix,ix])
  for(test in c('Mantel',if(j>1) 'Partial Mantel')) {
   fit <- if(test=='Mantel') mantel(x,y,method='pearson',permutations=perm) else
       mantel.partial(x,y,geo,method='pearson',permutations=perm)
   # Independently verify the reported correlation, using the same pair entries.
   rxy <- cor(as.vector(x),as.vector(y))
   expected <- if(test=='Mantel') rxy else {
    rxz <- cor(as.vector(x),as.vector(geo)); ryz <- cor(as.vector(y),as.vector(geo))
    (rxy-rxz*ryz)/sqrt((1-rxz^2)*(1-ryz^2))
   }
   stopifnot(abs(expected-unname(fit$statistic))<1e-12)
   k <- k+1L
   all_res[[k]] <- data.frame(test=test,variable=names(distances)[j],group=g,
    n_populations=length(ix),n_pairs=choose(length(ix),2),r=unname(fit$statistic),
    p=fit$signif,permutations=nrow(perm),seed=seed+gi-1L,
    control=if(test=='Mantel') 'None' else 'Geographic distance')
  }
 }
 message('Completed ',g,': ',length(ix),' populations')
}
res <- do.call(rbind,all_res)
# A single explicit family across the 52 displayed tests, for optional review.
res$p_BH_52 <- p.adjust(res$p,method='BH')
res$significance <- ifelse(res$p<0.001,'***',ifelse(res$p<0.01,'**',ifelse(res$p<0.05,'*','ns')))
write.csv(res,file.path(out,'mantel_results_long.csv'),row.names=FALSE,quote=TRUE)
saveRDS(perm_sets,file.path(out,'permutation_sets.rds'))
for(metric in c('r','p','p_BH_52')) {
 wide <- data.frame(variable=names(distances),check.names=FALSE)
 for(test in c('Mantel','Partial Mantel')) for(g in names(groups)) {
  rr <- res[res$test==test & res$group==g,]
  wide[[paste(test,g,sep='_')]] <- rr[[metric]][match(wide$variable,rr$variable)]
 }
 write.csv(wide,file.path(out,paste0('table_',metric,'.csv')),row.names=FALSE,na='')
}
capture.output(sessionInfo(),file=file.path(out,'sessionInfo.txt'))
write.table(data.frame(file=names(tools::md5sum(c(fst_path,file.path(base,files),
  file.path(base,c('pop_coord.txt','coords_climate_all.csv'))))),
  md5=unname(tools::md5sum(c(fst_path,file.path(base,files),
  file.path(base,c('pop_coord.txt','coords_climate_all.csv')))))),
  file.path(out,'input_md5.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
print(res[res$variable=='Geographic distance',c('group','n_populations','r','p')],row.names=FALSE)
