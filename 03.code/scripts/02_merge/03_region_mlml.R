#!/usr/bin/env Rscript
# name£º03_region_mlml.R
# input£º01_site/mC/hmC/C_matrix.tsv, annotion BED
# output£º{region}_mC.tsv, {region}_hmC.tsv, {region}_C.tsv
# usage£º03_region_mlml.R <site_dir> <anno_bed> <outdir> <region_name> <min_cpg> <agg>
# author:zhangyuxin,20260930

suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 6) {
  stop("Usage: mlml_03_region_level.R <site_dir> <anno_bed> <outdir> <region_name> <min_cpg> <agg>")
}

SITE_DIR    <- args[1]
ANNO_FILE   <- args[2]
OUTDIR      <- args[3]
REGION_NAME <- args[4]
MIN_CPG     <- as.integer(args[5])
AGG_METHOD  <- args[6]

dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

cat("Reading matrices ...\n")
mc  <- fread(file.path(SITE_DIR, "mC_matrix.tsv"))
hmc <- fread(file.path(SITE_DIR, "hmC_matrix.tsv"))
c_  <- fread(file.path(SITE_DIR, "C_matrix.tsv"))

anno <- fread(ANNO_FILE, header = FALSE,
              col.names = c("chr", "start", "end", "region_id"))

samples <- colnames(mc)[-(1:2)]
cat("Samples:", length(samples), "\n")
cat("Sites:", nrow(mc), "\n")
cat("Regions:", nrow(anno), "\n")

cpg_gr  <- GRanges(mc$chr, IRanges(mc$pos, mc$pos))
anno_gr <- GRanges(anno$chr, IRanges(anno$start, anno$end),
                   region_id = anno$region_id)

cat("Finding overlaps ...\n")
hits <- findOverlaps(anno_gr, cpg_gr)

agg_fun <- switch(AGG_METHOD,
  mean   = function(x) mean(x,   na.rm = TRUE),
  median = function(x) median(x, na.rm = TRUE),  
  function(x) mean(x, na.rm = TRUE)
)

aggregate_one <- function(mat) {
  res <- data.table(
    region_id = anno$region_id[hits@from],
    chr       = as.character(seqnames(anno_gr))[hits@from],
    start     = start(anno_gr)[hits@from],
    end       = end(anno_gr)[hits@from]
  )
  for (s in samples) {
    res[[s]] <- mat[[s]][hits@to]
  }
  agg <- res[, c(list(n_cpg = .N), lapply(.SD, agg_fun)),
             by = .(region_id, chr, start, end), .SDcols = samples]
  agg[n_cpg >= MIN_CPG]
}

cat("Aggregating mC ...\n")
mc_agg  <- aggregate_one(mc)
cat("Aggregating hmC ...\n")
hmc_agg <- aggregate_one(hmc)
cat("Aggregating C ...\n")
c_agg   <- aggregate_one(c_)

f_mc  <- file.path(OUTDIR, paste0(REGION_NAME, "_mC.tsv"))
f_hmc <- file.path(OUTDIR, paste0(REGION_NAME, "_hmC.tsv"))
f_c   <- file.path(OUTDIR, paste0(REGION_NAME, "_C.tsv"))

fwrite(mc_agg,  f_mc,  sep = "\t")
fwrite(hmc_agg, f_hmc, sep = "\t")
fwrite(c_agg,   f_c,   sep = "\t")

cat("Regions kept (mC):",  nrow(mc_agg),  "\n")
cat("Regions kept (hmC):", nrow(hmc_agg), "\n")
cat("Regions kept (C):",   nrow(c_agg),   "\n")
cat("Outputs:\n")
cat("  ", f_mc,  "\n")
cat("  ", f_hmc, "\n")
cat("  ", f_c,   "\n")