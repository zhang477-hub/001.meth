#!/usr/bin/env Rscript
# name£º02_window_mlml.R
# input£º01_site/mC_matrix.tsv, hmC_matrix.tsv, C_matrix.tsv
# output£ºwindow_{size}bp_mC.tsv, _hmC.tsv, _C.tsv,_nCpG.tsv
# usage£º02_window_mlml.R <site_dir> <outdir> <window> <step> <min_cpg> <min_cpg_per_sample>
# author:zhangyuxin,20260930
#modified:zhangyuxin,20261002

suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 6) {
  stop("Usage: 02_window_mlml.R <site_dir> <outdir> <window> <step> <min_cpg> <min_cpg_per_sample>")
}

SITE_DIR           <- args[1]
OUTDIR             <- args[2]
WINDOW             <- as.integer(args[3])
STEP               <- as.integer(args[4])
MIN_CPG            <- as.integer(args[5])
MIN_CPG_PER_SAMPLE <- as.integer(args[6])

dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

cat("==================================================\n")
cat("Start 02_window_mlml:", format(Sys.time()), "\n")
cat("SITE_DIR:            ", SITE_DIR, "\n")
cat("OUTDIR:              ", OUTDIR, "\n")
cat("WINDOW:              ", WINDOW, "\n")
cat("STEP:                ", STEP, "\n")
cat("MIN_CPG:             ", MIN_CPG, "\n")
cat("MIN_CPG_PER_SAMPLE:  ", MIN_CPG_PER_SAMPLE, "\n")
cat("==================================================\n")

cat("Reading matrices ...\n")
mc  <- fread(file.path(SITE_DIR, "mC_matrix.tsv"))
hmc <- fread(file.path(SITE_DIR, "hmC_matrix.tsv"))
c_  <- fread(file.path(SITE_DIR, "C_matrix.tsv"))

samples <- colnames(mc)[-(1:2)]
cat("Samples:", length(samples), "\n")
cat("Sites:", nrow(mc), "\n")

setkey(mc,  chr, pos)
setkey(hmc, chr, pos)
setkey(c_,  chr, pos)

chroms <- unique(mc$chr)
cat("Chromosomes:", length(chroms), "\n")

f_mc   <- file(file.path(OUTDIR, paste0("window_", WINDOW, "bp_mC.tsv")),  "w")
f_hmc  <- file(file.path(OUTDIR, paste0("window_", WINDOW, "bp_hmC.tsv")), "w")
f_c    <- file(file.path(OUTDIR, paste0("window_", WINDOW, "bp_C.tsv")),   "w")
f_ncpg <- file(file.path(OUTDIR, paste0("window_", WINDOW, "bp_nCpG.tsv")), "w")

header_meth <- paste(c("chr", "start", "end", "n_cpg", samples), collapse = "\t")
header_ncpg <- paste(c("chr", "start", "end", samples), collapse = "\t")

writeLines(header_meth, f_mc)
writeLines(header_meth, f_hmc)
writeLines(header_meth, f_c)
writeLines(header_ncpg, f_ncpg)

total_win <- 0

for (chr in chroms) {
  cat("  ", chr, "...\n")

  mc_chr  <- mc[chr == chr]
  hmc_chr <- hmc[chr == chr]
  c_chr   <- c_[chr == chr]

  if (nrow(mc_chr) == 0) next

  chr_len <- max(mc_chr$pos)
  starts <- seq(1, chr_len, by = STEP)

  for (st in starts) {
    en <- st + WINDOW - 1
    idx <- which(mc_chr$pos >= st & mc_chr$pos <= en)

    n_cpg_total <- length(idx)
    if (n_cpg_total < MIN_CPG) next

    mc_win  <- mc_chr[idx]
    hmc_win <- hmc_chr[idx]
    c_win   <- c_chr[idx]

    # non NA
    n_per_sample <- sapply(samples, function(s) {
      sum(!is.na(mc_win[[s]]) & !is.na(hmc_win[[s]]) & !is.na(c_win[[s]]))
    })

    # skip NA_CPG
    if (max(n_per_sample) < MIN_CPG_PER_SAMPLE) next

    row_common <- c(chr, st, en, n_cpg_total)
    row_mc  <- row_common
    row_hmc <- row_common
    row_c   <- row_common

    for (s in samples) {
      v1 <- mean(mc_win[[s]],  na.rm = TRUE)
      v2 <- mean(hmc_win[[s]], na.rm = TRUE)
      v3 <- mean(c_win[[s]],   na.rm = TRUE)

      row_mc  <- c(row_mc,  ifelse(is.nan(v1), "NA", sprintf("%.4f", v1)))
      row_hmc <- c(row_hmc, ifelse(is.nan(v2), "NA", sprintf("%.4f", v2)))
      row_c   <- c(row_c,   ifelse(is.nan(v3), "NA", sprintf("%.4f", v3)))
    }

    writeLines(paste(row_mc,  collapse = "\t"), f_mc)
    writeLines(paste(row_hmc, collapse = "\t"), f_hmc)
    writeLines(paste(row_c,   collapse = "\t"), f_c)
    writeLines(paste(c(chr, st, en, n_per_sample), collapse = "\t"), f_ncpg)

    total_win <- total_win + 1
  }
}

close(f_mc); close(f_hmc); close(f_c); close(f_ncpg)

cat("Windows written:", total_win, "\n")
cat("Outputs:\n")
cat("  ", file.path(OUTDIR, paste0("window_", WINDOW, "bp_mC.tsv")),  "\n")
cat("  ", file.path(OUTDIR, paste0("window_", WINDOW, "bp_hmC.tsv")), "\n")
cat("  ", file.path(OUTDIR, paste0("window_", WINDOW, "bp_C.tsv")),   "\n")
cat("  ", file.path(OUTDIR, paste0("window_", WINDOW, "bp_nCpG.tsv")), "\n")
cat("End 02_window_mlml:", format(Sys.time()), "\n")
