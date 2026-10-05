#!/usr/bin/env Rscript
# Per-chromosome RDA; no LD pruning. See run_RDA_campy.README.md.
main <- function(args) {
  usage <- paste(
    'Usage: Rscript run_RDA_campy.R input.raw environment.tsv [options]',
    '--pcs 0|2|3         Conditioning PCs (default 3: final conditioned model)',
    '--axes 6           Maximum constrained axes to scan',
    '--threads 1        data.table I/O threads, not RDA threads',
    '--batch-size 10000 Candidate SNPs per correlation batch',
    '--plot-snps 5000   Maximum SNP points in PDF (0 hides them)',
    '--id-column IID    Sample ID column: FID or IID',
    '--out-prefix PATH  Default: input path; append .pcN when N > 0', sep='\n')
  if (length(args) == 1L && args[1] %in% c('--help', '-h')) {
    cat(usage, '\n'); return(invisible(NULL))
  }
  if (length(args) < 2L) stop(usage, call.=FALSE)
  opt <- list(pcs='3', axes='6', threads='1', `batch-size`='10000',
              `plot-snps`='5000', `id-column`='IID', `out-prefix`=NULL)
  extra <- args[-c(1L, 2L)]
  if (length(extra) %% 2L) stop('Each option requires a value.\n', usage)
  if (length(extra)) for (i in seq.int(1L, length(extra), by=2L)) {
    key <- sub('^--', '', extra[i])
    if (!startsWith(extra[i], '--') || !key %in% names(opt)) stop('Unknown option: ', extra[i])
    opt[[key]] <- extra[i+1L]
  }
  for (key in c('pcs', 'axes', 'threads', 'batch-size', 'plot-snps')) {
    value <- suppressWarnings(as.numeric(opt[[key]]))
    minimum <- if (key %in% c('pcs', 'plot-snps')) 0 else 1
    if (length(value) != 1L || !is.finite(value) || value < minimum ||
        value > .Machine$integer.max || value != floor(value)) stop('Invalid --', key)
    opt[[key]] <- as.integer(value)
  }
  if (!opt$pcs %in% c(0L, 2L, 3L)) stop('--pcs must be 0, 2, or 3')
  if (!opt$`id-column` %in% c('FID', 'IID')) stop('--id-column must be FID or IID')
  if (!all(file.exists(args[1:2]))) stop('Input file does not exist')
  if (!requireNamespace('vegan', quietly=TRUE)) stop('Please install the vegan package')
  fast_io <- requireNamespace('data.table', quietly=TRUE)
  if (fast_io) data.table::setDTthreads(opt$threads)
  if (!fast_io) warning('data.table unavailable: using slower base R input/output')
  prefix <- opt$`out-prefix`
  if (is.null(prefix)) prefix <- paste0(args[1], if (opt$pcs) paste0('.pc', opt$pcs) else '')
  if (!dir.exists(dirname(prefix))) stop('Output directory does not exist: ', dirname(prefix))
  start <- proc.time()[[3L]]
  log_step <- function(...) cat(sprintf('[%.1f s] ', proc.time()[[3L]]-start), ..., '\n', sep='')
  predictors <- c('wc2.1_2.5m_srad_07', 'bio10', 'bio7', 'bio15', 'bio13', 'bio3')
  pcs <- if (opt$pcs) paste0('PC', seq_len(opt$pcs)) else character()
  log_step('Reading genotype file')
  if (fast_io) {
    raw <- data.table::fread(args[1], header=TRUE, na.strings=c('NA', 'NaN'),
                            colClasses=list(character=c('FID', 'IID')), check.names=FALSE)
  } else {
    raw <- read.table(args[1], header=TRUE, sep='', check.names=FALSE,
                      colClasses=c(FID='character', IID='character'), comment.char='',
                      quote='', stringsAsFactors=FALSE)
  }
  if (ncol(raw) < 7L || !identical(names(raw)[1:6],
      c('FID', 'IID', 'PAT', 'MAT', 'SEX', 'PHENOTYPE'))) stop('Expected standard PLINK additive .raw header')
  if (anyDuplicated(names(raw))) stop('Duplicate column/SNP names')
  keep <- which(!raw[['IID']] %in% c('xubo1407', 'XB_DR_C'))
  ids <- raw[[opt$`id-column`]][keep]
  if (anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids)) stop('Invalid/duplicate sample IDs; consider --id-column IID')
  if (length(ids) < 3L) stop('Too few samples')
  ordering <- order(ids, decreasing=TRUE)
  ids <- ids[ordering]
  keep <- keep[ordering]
  raw <- if (fast_io) raw[keep] else raw[keep, , drop=FALSE]
  env <- read.delim(args[2], header=TRUE, check.names=FALSE,
                    colClasses=c(sample='character'), stringsAsFactors=FALSE)
  required <- c('sample', predictors, pcs)
  missing <- setdiff(required, names(env))
  if (length(missing)) stop('Missing environmental columns: ', paste(missing, collapse=', '))
  if (anyNA(env$sample) || anyDuplicated(env$sample)) stop('Invalid/duplicate environment sample IDs')
  idx <- match(ids, env$sample)
  if (anyNA(idx)) stop('No environmental row for: ', paste(ids[is.na(idx)], collapse=', '))
  env <- env[idx, , drop=FALSE]
  stopifnot(identical(ids, env$sample))
  variables <- c(predictors, pcs)
  valid <- vapply(env[variables], function(x) {
    is.numeric(x) && all(is.finite(x)) && sd(x) > 0
  }, logical(1))
  if (!all(valid)) stop('Non-numeric, missing, or constant predictors: ', paste(variables[!valid], collapse=', '))
  pred <- as.data.frame(scale(env[variables]))
  rownames(pred) <- ids
  design <- cbind(Intercept=1, as.matrix(pred))
  if (qr(design)$rank < ncol(design)) stop('Predictors/PCs are linearly dependent; review the model')
  if (nrow(design) <= ncol(design)) stop('No residual degrees of freedom')
  log_step('Checking SNPs and imputing missing genotypes by mode')
  snp_columns <- seq.int(7L, ncol(raw))
  retained <- logical(length(snp_columns))
  missing_count <- 0
  all_missing <- constant <- 0L
  for (k in seq_along(snp_columns)) {
    j <- snp_columns[k]
    x <- raw[[j]]
    if (!is.numeric(x) && !(is.logical(x) && all(is.na(x)))) stop('Non-numeric genotypes: ', names(raw)[j])
    nas <- which(is.na(x))
    missing_count <- missing_count + length(nas)
    if (length(nas) == length(x)) { all_missing <- all_missing+1L; next }
    if (any(!is.na(x) & !(x %in% c(0, 1, 2)))) stop('Expected hard-call 0/1/2/NA genotypes: ', names(raw)[j])
    counts <- tabulate(as.integer(x)+1L, nbins=3L)
    if (sum(counts > 0L) < 2L) { constant <- constant+1L; next }
    retained[k] <- TRUE
    if (length(nas)) {
      # Same tie-breaking as the original table/which.max: smallest dosage wins.
      mode <- which.max(counts)-1L
      if (fast_io) data.table::set(raw, i=nas, j=j, value=mode) else {
        x[nas] <- mode; raw[[j]] <- x
      }
    }
  }
  selected <- snp_columns[retained]
  if (length(selected) < 2L) stop('Too few polymorphic SNPs')
  gen <- if (fast_io) as.matrix(raw[, selected, with=FALSE]) else as.matrix(raw[, selected, drop=FALSE])
  rownames(gen) <- ids
  rm(raw)
  invisible(gc())
  log_step('Samples: ', nrow(gen), '; SNPs: ', ncol(gen), '; original missing calls: ', missing_count,
           '; removed all-missing: ', all_missing, '; removed constant: ', constant)
  rhs <- paste(predictors, collapse=' + ')
  if (length(pcs)) rhs <- paste(rhs, '+ Condition(', paste(pcs, collapse=' + '), ')')
  formula <- as.formula(paste('gen ~', rhs))
  log_step('Fitting RDA with ', opt$pcs, ' conditioning PCs')
  fit <- vegan::rda(formula, data=pred, scale=TRUE)
  rank <- fit$CCA$rank
  if (is.null(rank) || rank < 1L) stop('No constrained axes')
  # Use every constrained eigenvalue, independently of --axes (scan limit).
  eigenvalues <- fit$CCA$eig
  constrained_inertia <- sum(eigenvalues)
  axis_variance <- data.frame(
    axis=names(eigenvalues),
    eigenvalue=unname(eigenvalues),
    percent_total=100 * unname(eigenvalues) / fit$tot.chi,
    percent_constrained=100 * unname(eigenvalues) / constrained_inertia,
    row.names=NULL)
  write.table(axis_variance, paste0(prefix, '_rda_axis_variance.tsv'),
              sep='\t', quote=FALSE, row.names=FALSE)
  log_step('RDA axis explained variance (percent):')
  print(head(axis_variance, 2L), row.names=FALSE, digits=6)
  axes <- seq_len(min(opt$axes, rank))
  if (length(axes) < opt$axes) warning('Requested axes exceed constrained rank; using ', length(axes))
  loadings <- as.matrix(vegan::scores(fit, choices=axes, display='species', scaling=2))
  sites <- as.matrix(vegan::scores(fit, choices=axes, display='sites', scaling=2))
  write.table(sites, paste0(prefix, 'rda_ind.txt'), quote=FALSE, col.names=NA)
  candidates <- function(z) {
    parts <- lapply(seq_along(axes), function(k) {
      x <- loadings[, k]
      selected <- which(abs(x-mean(x)) > z*sd(x))
      data.frame(axis=rep(axes[k], length(selected)), snp=rownames(loadings)[selected],
                 loading=unname(x[selected]), stringsAsFactors=FALSE)
    })
    result <- do.call(rbind, parts)
    # Preserve original first-axis deduplication.
    result[!duplicated(result$snp), , drop=FALSE]
  }
  thresholds <- c(3.5, 4.0)
  candidate_sets <- lapply(thresholds, candidates)
  unique_snps <- unique(unlist(lapply(candidate_sets, function(x) x$snp), use.names=FALSE))
  correlations <- matrix(NA_real_, length(unique_snps), length(predictors),
                         dimnames=list(unique_snps, predictors))
  environmental_matrix <- as.matrix(pred[predictors])
  log_step('Computing correlations for ', length(unique_snps), ' unique candidates')
  if (length(unique_snps)) for (first in seq.int(1L, length(unique_snps), by=opt$`batch-size`)) {
    rows <- seq.int(first, min(first+opt$`batch-size`-1, length(unique_snps)))
    correlations[rows, ] <- cor(gen[, unique_snps[rows], drop=FALSE], environmental_matrix)
  }
  for (k in seq_along(thresholds)) {
    cand <- candidate_sets[[k]]
    cc <- correlations[match(cand$snp, unique_snps), , drop=FALSE]
    best <- if (nrow(cand)) max.col(abs(cc), ties.method='first') else integer()
    rownames(cand) <- NULL
    rownames(cc) <- NULL
    cand <- cbind(cand, as.data.frame(cc))
    cand$predictor <- predictors[best]
    cand$correlation <- if (nrow(cand)) abs(cc[cbind(seq_len(nrow(cand)), best)]) else numeric()
    path <- paste0(prefix, 'rdaz', sprintf('%.1f', thresholds[k]), '.txt')
    # Keep the legacy row-number column for downstream compatibility.
    if (fast_io) data.table::fwrite(cand, path, sep='\t', quote=FALSE, row.names=TRUE) else {
      write.table(cand, path, quote=FALSE, sep='\t', row.names=TRUE)
    }
    log_step('Threshold ', thresholds[k], ': ', nrow(cand), ' unique candidate SNPs')
  }
  if (rank >= 2L) {
    log_step('Writing plot; maximum SNP points: ', opt$`plot-snps`)
    local({
      pdf(paste0(prefix, '_rda_plot.pdf'))
      on.exit(dev.off())
      plot(fit, type='n', display=c('sites', 'bp'), choices=1:2, scaling=3)
      if (opt$`plot-snps` > 0L) {
        sp <- as.matrix(vegan::scores(fit, choices=1:2, display='species', scaling=3))
        set.seed(1)
        pick <- if (nrow(sp) > opt$`plot-snps`) sample.int(nrow(sp), opt$`plot-snps`) else seq_len(nrow(sp))
        points(sp[pick, , drop=FALSE], pch=20, cex=0.5, col='gray70')
      }
      group <- if ('group' %in% names(env)) factor(env$group) else factor(rep('Samples', nrow(env)))
      palette <- setNames(grDevices::hcl.colors(nlevels(group), 'Dark 3'), levels(group))
      known <- c(CPO='orange', CYU='steelblue', CMA='darkgreen', N='grey')
      common <- intersect(names(palette), names(known))
      palette[common] <- known[common]
      points(fit, choices=1:2, display='sites', pch=21, cex=1.3, col='gray32',
             scaling=3, bg=unname(palette[as.character(group)]))
      text(fit, choices=1:2, scaling=3, display='bp', col='darkred', cex=0.8)
      legend('bottomright', legend=levels(group), bty='n', col='gray32', pch=21,
             pt.bg=unname(palette), cex=0.8)
    })
  } else warning('Only one constrained axis: two-dimensional PDF skipped')
  report <- capture.output({
    cat('Command arguments:\n'); print(args)
    cat('\nModel:\n'); print(formula)
    cat('\nSamples:', nrow(gen), ' SNPs:', ncol(gen), '\n')
    cat('Removed all-missing:', all_missing, ' Constant:', constant, '\n')
    cat('Scanned constrained axes:', axes, '\n')
    cat('Axes are selected by count, not by a significance test.\n')
    cat('Candidate thresholds are loading SD cutoffs, not corrected P-values.\n')
    cat('Candidate correlations are marginal Pearson correlations, not PC-adjusted effects.\n')
    cat('\nR-squared:\n'); print(vegan::RsquareAdj(fit))
    cat('\nVIF:\n'); print(vegan::vif.cca(fit))
    cat('\nConstrained eigenvalues:\n'); print(fit$CCA$eig)
    cat('\nAxis explained variance (percent):\n')
    cat('percent_total: eigenvalue / total SNP inertia, before PC conditioning.\n')
    cat('percent_constrained: eigenvalue / sum of ALL constrained eigenvalues.\n')
    cat('With conditioning PCs, constrained inertia is the predictor contribution after PC adjustment.\n')
    cat('Total inertia:', format(fit$tot.chi, digits=12), '\n')
    cat('Constrained inertia:', format(constrained_inertia, digits=12), '\n')
    print(axis_variance, row.names=FALSE, digits=6)
    cat('First', min(2L, nrow(axis_variance)), 'constrained axes combined (%):\n')
    print(colSums(head(axis_variance[c('percent_total', 'percent_constrained')], 2L)), digits=6)
    cat('\nElapsed seconds:', proc.time()[[3L]]-start, '\n')
    cat('\nSession:\n'); print(sessionInfo())
  })
  writeLines(report, paste0(prefix, '_rda_summary.txt'))
  log_step('Finished; output prefix: ', prefix)
}
if (sys.nframe() == 0L) main(commandArgs(trailingOnly=TRUE))
