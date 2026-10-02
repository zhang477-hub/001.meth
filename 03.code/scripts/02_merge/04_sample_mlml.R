#!/usr/bin/env Rscript
# name£º04_sample_mlml.R
# input£º01_site/mC/hmC/C_matrix.tsv
# output£ºglobal_meth.tsv
# usage£º04_sample_mlml.R <site_dir> <outfile>
# author:zhangyuxin,20260930

suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: mlml_04_sample_level.R <site_dir> <outfile>")
}

SITE_DIR <- args[1]
OUTFILE  <- args[2]

dir.create(dirname(OUTFILE), recursive = TRUE, showWarnings = FALSE)

cat("Reading matrices ...\n")
mc  <- fread(file.path(SITE_DIR, "mC_matrix.tsv"))
hmc <- fread(file.path(SITE_DIR, "hmC_matrix.tsv"))
c_  <- fread(file.path(SITE_DIR, "C_matrix.tsv"))

samples <- colnames(mc)[-(1:2)]
cat("Samples:", length(samples), "\n")

res <- data.table(
  sample     = samples,
  global_mC  = sapply(samples, function(s) mean(mc[[s]],  na.rm = TRUE)),
  global_hmC = sapply(samples, function(s) mean(hmc[[s]], na.rm = TRUE)),
  global_C   = sapply(samples, function(s) mean(c_[[s]],  na.rm = TRUE)),
  n_sites    = sapply(samples, function(s) sum(!is.na(mc[[s]])))
)

fwrite(res, OUTFILE, sep = "\t")

cat("Output:", OUTFILE, "\n")
cat("Samples:", nrow(res), "\n")