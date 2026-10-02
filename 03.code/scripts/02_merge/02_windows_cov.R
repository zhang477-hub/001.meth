#!/usr/bin/env Rscript
# name£º02_windows_cov.R
# input£º01_site/beta_matrix.tsv, 01_site/coverage_matrix.tsv
# output£ºwindow_{size}bp.tsv,window_{size}bp_nCpG.tsv
# usage£º02_windows_cov.R <beta_matrix> <cov_matrix> <outfile> <window> <step> <min_cpg> <min_cov> <min_cpg_per_sample>
# author:zhangyuxin,20260929
#modified:zhangyuxin,20261002,add n_cpg

suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 8) {
  stop("Usage: 02_windows_cov.R <beta> <cov> <outdir> <window> <step> <min_cpg> <min_cov> <min_cpg_per_sample>")
}

BETA_FILE          <- args[1]
COV_FILE           <- args[2]
OUTDIR             <- args[3]
WINDOW             <- as.integer(args[4])
STEP               <- as.integer(args[5])
MIN_CPG            <- as.integer(args[6])
MIN_COV            <- as.integer(args[7])
MIN_CPG_PER_SAMPLE <- as.integer(args[8])

dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

cat("==================================================\n")
cat("Start 02_windows_cov:", format(Sys.time()), "\n")
cat("BETA_FILE:           ", BETA_FILE, "\n")
cat("COV_FILE:            ", COV_FILE, "\n")
cat("OUTDIR:              ", OUTDIR, "\n")
cat("WINDOW:              ", WINDOW, "\n")
cat("STEP:                ", STEP, "\n")
cat("MIN_CPG:             ", MIN_CPG, "\n")
cat("MIN_COV:             ", MIN_COV, "\n")
cat("MIN_CPG_PER_SAMPLE:  ", MIN_CPG_PER_SAMPLE, "\n")
cat("==================================================\n")

cat("Reading beta ...\n")
beta <- fread(BETA_FILE)
cat("Reading cov ...\n")
cov  <- fread(COV_FILE)

samples <- colnames(beta)[-(1:2)]
cat("Samples:", length(samples), "\n")
cat("Sites:", nrow(beta), "\n")

setkey(beta, chr, pos)
setkey(cov, chr, pos)

chroms <- unique(beta$chr)
cat("Chromosomes:", length(chroms), "\n")

# output
OUTFILE <- file.path(OUTDIR, paste0("window_", WINDOW, "bp.tsv"))
NCPG_FILE <- file.path(OUTDIR, paste0("window_", WINDOW, "bp_nCpG.tsv"))

fout  <- file(OUTFILE, "w")
f_ncpg <- file(NCPG_FILE, "w")

writeLines(paste(c("chr", "start", "end", "n_cpg", "mean_cov", samples),
                 collapse = "\t"), fout)
writeLines(paste(c("chr", "start", "end", samples),
                 collapse = "\t"), f_ncpg)

total_win <- 0

for (chr in chroms) {
  cat("  ", chr, "...\n")

  b_chr <- beta[chr == chr]
  c_chr <- cov[chr == chr]

  if (nrow(b_chr) == 0) next

  chr_len <- max(b_chr$pos)
  starts <- seq(1, chr_len, by = STEP)

  for (st in starts) {
    en <- st + WINDOW - 1
    idx <- which(b_chr$pos >= st & b_chr$pos <= en)

    n_cpg_total <- length(idx)
    if (n_cpg_total < MIN_CPG) next

    b_win <- b_chr[idx]
    c_win <- c_chr[idx]

    n_per_sample <- sapply(samples, function(s) {
      m  <- b_win[[s]]
      cv <- c_win[[s]]
      sum(!is.na(m) & !is.na(cv) & cv >= MIN_COV)
    })

    if (max(n_per_sample) < MIN_CPG_PER_SAMPLE) next

    row <- list(chr = chr, start = st, end = en, n_cpg = n_cpg_total)

    cov_vals <- c()
    for (s in samples) {
      m  <- b_win[[s]]
      cv <- c_win[[s]]
      valid <- !is.na(m) & !is.na(cv) & cv >= MIN_COV
      if (sum(valid) == 0) {
        row[[s]] <- NA_real_
      } else {
        row[[s]] <- sum(m[valid] * cv[valid]) / sum(cv[valid])
      }
      cov_vals <- c(cov_vals, mean(cv, na.rm = TRUE))
    }
    row$mean_cov <- mean(cov_vals, na.rm = TRUE)

    line <- paste(c(row$chr, row$start, row$end, row$n_cpg,
                    sprintf("%.2f", row$mean_cov),
                    sapply(samples, function(s) {
                      v <- row[[s]]
                      if (is.na(v)) "NA" else sprintf("%.4f", v)
                    })), collapse = "\t")
    writeLines(line, fout)
    writeLines(paste(c(chr, st, en, n_per_sample), collapse = "\t"), f_ncpg)

    total_win <- total_win + 1
  }
}

close(fout)
close(f_ncpg)

cat("==================================================\n")
cat("Windows written:", total_win, "\n")
cat("Outputs:\n")
cat("  ", OUTFILE,   "\n")
cat("  ", NCPG_FILE, "\n")
cat("End 02_windows_cov:", format(Sys.time()), "\n")
cat("==================================================\n")