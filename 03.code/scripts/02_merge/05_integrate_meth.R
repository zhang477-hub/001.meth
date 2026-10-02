#!/usr/bin/env Rscript
# name：05_integrate_meth.R
# input：02_window/window_{size}bp_{mC,hmC,C}.tsv, nCpG.tsv, 04_sample/global_meth.tsv, pheno.tsv
# output：meth_1kb_impute.rda
# usage：05_integrate_meth.R <window_dir> <sample_qc> <pheno> <outdir> <window_size> [options]
# author:zhangyuxin,20261002

suppressPackageStartupMessages({
  library(data.table)
})

# ============================================================
# 参数解析
# ============================================================
parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 5) {
    cat("Usage: 05_integrate_meth.R <window_dir> <sample_qc> <pheno> <outdir> <window_size> [options]\n")
    cat("  window_dir   : 02_window/window_{W}bp_mC.tsv 等\n")
    cat("  sample_qc    : 04_sample/global_meth.tsv\n")
    cat("  pheno        :  sample_id, group, batch \n")
    cat("  outdir       : output\n")
    cat("  window_size  : eg 1000\n")
    cat("Options:\n")
    cat("  --min-conv 0.95        coversion（ 0.95）\n")
    cat("  --min-cov 5            ave coverage（ 5）\n")
    cat("  --min-global-mC 0.3    global meth min（ 0.3）\n")
    cat("  --min-global-hmC 0.15  global hmeth max（ 0.15）\n")
    cat("  --min-cpg 5            min CpG number（ 5）\n")
    cat("  --max-missing 0.5      NA（ 0.5）\n")
    cat("  --impute knn           input method：knn / mean / median / missforest / none（ knn）\n")
    cat("  --normalize quantile   sclaed：quantile / combat / none（ none）\n")
    cat("  --batch-col batch      （ batch）\n")
    cat("  --group-col group      （ group）\n")
    quit(status = 0)
  }
  opts <- list(
    window_dir      = args[1],
    sample_qc       = args[2],
    pheno           = args[3],
    outdir          = args[4],
    window_size     = as.integer(args[5]),
    min_conv        = 0.95,
    min_cov         = 5,
    min_global_mC   = 0.3,
    min_global_hmC  = 0.15,
    min_cpg         = 5,
    max_missing     = 0.5,
    impute          = "knn",
    normalize       = "none",
    batch_col       = "batch",
    group_col       = "group"
  )
  i <- 6
  while (i <= length(args)) {
    key <- args[i]
    val <- if (i < length(args)) args[i + 1] else NULL
    if (key == "--min-conv")         { opts$min_conv       <- as.numeric(val); i <- i + 2 }
    else if (key == "--min-cov")     { opts$min_cov        <- as.numeric(val); i <- i + 2 }
    else if (key == "--min-global-mC")  { opts$min_global_mC  <- as.numeric(val); i <- i + 2 }
    else if (key == "--min-global-hmC") { opts$min_global_hmC <- as.numeric(val); i <- i + 2 }
    else if (key == "--min-cpg")     { opts$min_cpg        <- as.integer(val); i <- i + 2 }
    else if (key == "--max-missing") { opts$max_missing    <- as.numeric(val); i <- i + 2 }
    else if (key == "--impute")      { opts$impute         <- val; i <- i + 2 }
    else if (key == "--normalize")   { opts$normalize      <- val; i <- i + 2 }
    else if (key == "--batch-col")   { opts$batch_col      <- val; i <- i + 2 }
    else if (key == "--group-col")   { opts$group_col      <- val; i <- i + 2 }
    else { stop(paste("Unknown option:", key)) }
  }
  opts
}

opts <- parse_args()

W        <- opts$window_size
WDIR     <- opts$window_dir
OUTDIR   <- opts$outdir
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

LOG_FILE <- file.path(OUTDIR, "05_integrate_meth.log")
log_con  <- file(LOG_FILE, open = "at")
sink(log_con, append = TRUE)
sink(log_con, append = TRUE, type = "message")

cat("==================================================\n")
cat("Start 05_integrate_meth:", format(Sys.time()), "\n")
cat("window_dir:  ", WDIR, "\n")
cat("window_size: ", W, "\n")
cat("outdir:      ", OUTDIR, "\n")
cat("impute:      ", opts$impute, "\n")
cat("normalize:   ", opts$normalize, "\n")
cat("==================================================\n")

#input

cat("\n[1] Reading matrices ...\n")

mc_file  <- file.path(WDIR, paste0("window_", W, "bp_mC.tsv"))
hmc_file <- file.path(WDIR, paste0("window_", W, "bp_hmC.tsv"))
c_file   <- file.path(WDIR, paste0("window_", W, "bp_C.tsv"))
ncpg_file <- file.path(WDIR, paste0("window_", W, "bp_nCpG.tsv"))

for (f in c(mc_file, hmc_file, c_file, ncpg_file)) {
  if (!file.exists(f)) stop(paste("File not found:", f))
}

mc   <- fread(mc_file)
hmc  <- fread(hmc_file)
c_   <- fread(c_file)
ncpg <- fread(ncpg_file)

cat("  mC  dims:", dim(mc),   "\n")
cat("  hmC dims:", dim(hmc),  "\n")
cat("  C   dims:", dim(c_),   "\n")
cat("  nCpG dims:", dim(ncpg), "\n")

# check indentity
if (!identical(mc$chr, hmc$chr) || !identical(mc$start, hmc$start)) {
  stop("mC and hmC windows are not aligned")
}

regions <- paste(mc$chr, mc$start, mc$end, sep = ":")
samples <- colnames(mc)[-(1:4)]  # 去掉 chr/start/end/n_cpg

cat("  Total regions:", length(regions), "\n")
cat("  Total samples:", length(samples), "\n")

# sample,phenotsv
cat("\n[2] Reading QC and phenotype ...\n")

qc    <- fread(opts$sample_qc)
pheno <- fread(opts$pheno)

cat("  QC samples:", nrow(qc), "\n")
cat("  Pheno samples:", nrow(pheno), "\n")

# processing QC
qc_cols <- colnames(qc)
cat("  QC columns:", paste(qc_cols, collapse = ", "), "\n")

# cov conversion
get_col <- function(df, candidates, default = NA) {
  for (c in candidates) if (c %in% colnames(df)) return(c)
  return(default)
}
conv_col <- get_col(qc, c("conv_eff", "conversion_efficiency", "conv_eff_pct"))
cov_col  <- get_col(qc, c("mean_cov", "avg_coverage", "coverage"))

# filter low_quality
cat("\n[3] Filtering low-quality samples ...\n")

keep_samples <- samples

# 3.1 QC 
qc_keep <- qc$sample
if (!is.na(conv_col)) {
  qc_keep <- intersect(qc_keep,
                       qc$sample[qc[[conv_col]] >= opts$min_conv * 100])
  cat("  After conv_eff >= ", opts$min_conv * 100, "%: ", length(qc_keep), "\n", sep = "")
}
if (!is.na(cov_col)) {
  qc_keep <- intersect(qc_keep,
                       qc$sample[qc[[cov_col]] >= opts$min_cov])
  cat("  After mean_cov >= ", opts$min_cov, ": ", length(qc_keep), "\n", sep = "")
}
keep_samples <- intersect(keep_samples, qc_keep)

# 3.2 meth
if ("global_mC" %in% colnames(qc)) {
  keep_samples <- intersect(keep_samples,
                            qc$sample[qc$global_mC >= opts$min_global_mC])
  cat("  After global_mC >= ", opts$min_global_mC, ": ", length(keep_samples), "\n", sep = "")
}
if ("global_hmC" %in% colnames(qc)) {
  keep_samples <- intersect(keep_samples,
                            qc$sample[qc$global_hmC <= opts$min_global_hmC])
  cat("  After global_hmC <= ", opts$min_global_hmC, ": ", length(keep_samples), "\n", sep = "")
}

# 3.3 pheno
keep_samples <- intersect(keep_samples, pheno$sample_id)

cat("  Final samples kept:", length(keep_samples), "/", length(samples), "\n")

if (length(keep_samples) < 3) {
  stop("Too few samples after filtering")
}

# t
cat("\n[4] Transposing to sample × region matrix ...\n")

mat_mc  <- as.matrix(mc[,  ..keep_samples])
mat_hmc <- as.matrix(hmc[, ..keep_samples])
mat_c   <- as.matrix(c_[,  ..keep_samples])
mat_n   <- as.matrix(ncpg[, ..keep_samples])

rownames(mat_mc)  <- regions
rownames(mat_hmc) <- regions
rownames(mat_c)   <- regions
rownames(mat_n)   <- regions
colnames(mat_mc)  <- keep_samples
colnames(mat_hmc) <- keep_samples
colnames(mat_c)   <- keep_samples
colnames(mat_n)   <- keep_samples

cat("  mC  dims:", dim(mat_mc), "\n")
cat("  hmC dims:", dim(mat_hmc), "\n")
cat("  C   dims:", dim(mat_c), "\n")

# filter qujian
cat("\n[5] Filtering low-quality regions ...\n")

# 5.1 CpG num filter
frac_ok_cpg <- rowMeans(mat_n >= opts$min_cpg, na.rm = TRUE)
keep_region_cpg <- which(frac_ok_cpg >= 0.8)
cat("  After CpG >= ", opts$min_cpg, " in >=80% samples: ", length(keep_region_cpg), "\n", sep = "")

# 5.2 NA filter
missing_rate <- rowMeans(is.na(mat_mc))
keep_region_miss <- which(missing_rate <= opts$max_missing)
cat("  After missing rate <= ", opts$max_missing, ": ", length(keep_region_miss), "\n", sep = "")

keep_regions <- intersect(keep_region_cpg, keep_region_miss)
cat("  Final regions kept:", length(keep_regions), "/", length(regions), "\n")

if (length(keep_regions) < 100) {
  stop("Too few regions after filtering")
}

# apply
mat_mc  <- mat_mc[keep_regions,  , drop = FALSE]
mat_hmc <- mat_hmc[keep_regions, , drop = FALSE]
mat_c   <- mat_c[keep_regions,   , drop = FALSE]
mat_n   <- mat_n[keep_regions,   , drop = FALSE]

# input complementry
cat("\n[6] Imputing missing values (", opts$impute, ") ...\n", sep = "")

impute_matrix <- function(mat, method) {
  n_na <- sum(is.na(mat))
  cat("  Missing values:", n_na, "/", length(mat),
      sprintf("(%.2f%%)", n_na / length(mat) * 100), "\n")

  if (n_na == 0) {
    cat("  No missing values, skip imputation\n")
    return(mat)
  }

  if (method == "none") {
    cat("  Imputation skipped\n")
    return(mat)
  }

  mat_t <- t(mat)

  if (method == "mean") {
    for (j in seq_len(ncol(mat_t))) {
      v <- mat_t[, j]
      if (any(is.na(v))) {
        mat_t[is.na(v), j] <- mean(v, na.rm = TRUE)
      }
    }
  } else if (method == "median") {
    for (j in seq_len(ncol(mat_t))) {
      v <- mat_t[, j]
      if (any(is.na(v))) {
        mat_t[is.na(v), j] <- median(v, na.rm = TRUE)
      }
    }
  } else if (method == "knn") {
    if (!requireNamespace("impute", quietly = TRUE)) {
      stop("Package 'impute' not installed. Run: BiocManager::install('impute')")
    }
    library(impute)
    mat_t <- impute.knn(mat_t, k = 10)$data
  } else if (method == "missforest") {
    if (!requireNamespace("missForest", quietly = TRUE)) {
      stop("Package 'missForest' not installed. Run: install.packages('missForest')")
    }
    library(missForest)
    set.seed(42)
    mat_t <- missForest(mat_t)$ximp
  } else {
    stop(paste("Unknown impute method:", method))
  }

  t(mat_t)
}

mat_mc  <- impute_matrix(mat_mc,  opts$impute)
mat_hmc <- impute_matrix(mat_hmc, opts$impute)
mat_c   <- impute_matrix(mat_c,   opts$impute)

# scaled
cat("\n[7] Normalizing (", opts$normalize, ") ...\n", sep = "")

normalize_matrix <- function(mat, method, pheno, batch_col) {
  if (method == "none") {
    cat("  Normalization skipped\n")
    return(mat)
  }

  if (method == "quantile") {
    if (!requireNamespace("preprocessCore", quietly = TRUE)) {
      stop("Package 'preprocessCore' not installed. Run: BiocManager::install('preprocessCore')")
    }
    library(preprocessCore)
    mat_norm <- normalize.quantiles(mat)
    dimnames(mat_norm) <- dimnames(mat)
    return(mat_norm)
  }

  if (method == "combat") {
    if (!requireNamespace("sva", quietly = TRUE)) {
      stop("Package 'sva' not installed. Run: BiocManager::install('sva')")
    }
    library(sva)
    if (!batch_col %in% colnames(pheno)) {
      stop(paste("Batch column not found:", batch_col))
    }
    batch <- pheno[[batch_col]][match(colnames(mat), pheno$sample_id)]
    if (any(is.na(batch))) {
      stop("Some samples have no batch info")
    }
    mod <- model.matrix(~ 1, data = pheno)
    mat_norm <- ComBat(mat, batch = batch, mod = mod, par.prior = TRUE)
    return(mat_norm)
  }

  stop(paste("Unknown normalize method:", method))
}

mat_mc  <- normalize_matrix(mat_mc,  opts$normalize, pheno, opts$batch_col)
mat_hmc <- normalize_matrix(mat_hmc, opts$normalize, pheno, opts$batch_col)
mat_c   <- normalize_matrix(mat_c,   opts$normalize, pheno, opts$batch_col)

# output
cat("\n[8] Preparing output object ...\n")

pheno_keep <- pheno[match(colnames(mat_mc), pheno$sample_id), ]
rownames(pheno_keep) <- pheno_keep$sample_id

region_info <- data.table(
  region_id = rownames(mat_mc),
  chr       = mc$chr[keep_regions],
  start     = mc$start[keep_regions],
  end       = mc$end[keep_regions]
)

meth_1kb <- list(
  mC        = mat_mc,         
  hmC       = mat_hmc,         
  C         = mat_c,           
  nCpG      = mat_n,           
  pheno     = pheno_keep,      
  region    = region_info,     
  params    = opts             
)

# downstream input
meth_1kb$mC_t   <- t(mat_mc)
meth_1kb$hmC_t  <- t(mat_hmc)
meth_1kb$C_t    <- t(mat_c)

#save
cat("\n[9] Saving .rda ...\n")

OUTFILE <- file.path(OUTDIR, "meth_1kb_impute.rda")
save(meth_1kb, file = OUTFILE, compress = "xz")

cat("  Saved:", OUTFILE, "\n")
cat("  File size:", round(file.size(OUTFILE) / 1024^2, 2), "MB\n")

# summary
cat("\n==================================================\n")
cat("Summary:\n")
cat("  Samples in:   ", length(samples), "\n")
cat("  Samples kept: ", length(keep_samples), "\n")
cat("  Regions in:   ", length(regions), "\n")
cat("  Regions kept: ", length(keep_regions), "\n")
cat("  Impute:       ", opts$impute, "\n")
cat("  Normalize:    ", opts$normalize, "\n")
cat("  Output:       ", OUTFILE, "\n")
cat("End 05_integrate_meth:", format(Sys.time()), "\n")
cat("==================================================\n")

sink(type = "message")
sink()
close(log_con)