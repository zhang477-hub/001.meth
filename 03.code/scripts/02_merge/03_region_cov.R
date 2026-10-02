#!/usr/bin/env Rscript
# name£º03_region_cov.R
# input£º01_site/beta_matrix.tsv, annotion BED
# output£º{region_type}_matrix.tsv
# usage£º03_region_cov.R <beta_matrix> <anno_bed> <outfile> <min_cpg> <agg_method>
# author:zhangyuxin,20260929

suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
  stop("Usage: 03_region_level.R <beta> <anno_bed> <outfile> <min_cpg> <agg>")
}

BETA_FILE  <- args[1]
ANNO_FILE  <- args[2]
OUTFILE    <- args[3]
MIN_CPG    <- as.integer(args[4])
AGG_METHOD <- args[5]   # mean / median / weighted

dir.create(dirname(OUTFILE), recursive = TRUE, showWarnings = FALSE)

cat("Reading beta ...\n")
beta <- fread(BETA_FILE)
cat("Reading annotation ...\n")
anno <- fread(ANNO_FILE, header = FALSE,
              col.names = c("chr", "start", "end", "region_id"))

samples <- colnames(beta)[-(1:2)]
cat("Samples:", length(samples), "\n")
cat("Sites:", nrow(beta), "\n")
cat("Regions:", nrow(anno), "\n")

# GRanges
cpg_gr  <- GRanges(beta$chr, IRanges(beta$pos, beta$pos))
anno_gr <- GRanges(anno$chr, IRanges(anno$start, anno$end),
                   region_id = anno$region_id)

cat("Finding overlaps ...\n")
hits <- findOverlaps(anno_gr, cpg_gr)

# results
res <- data.table(
  region_id = anno$region_id[hits@from],
  chr       = as.character(seqnames(anno_gr))[hits@from],
  start     = start(anno_gr)[hits@from],
  end       = end(anno_gr)[hits@from]
)

for (s in samples) {
  res[[s]] <- beta[[s]][hits@to]
}

# affregated
cat("Aggregating by", AGG_METHOD, "...\n")

agg_fun <- switch(AGG_METHOD,
  mean     = function(x) mean(x, na.rm = TRUE),
  median   = function(x) median(x, na.rm = TRUE),
  weighted = function(x) mean(x, na.rm = TRUE),
  function(x) mean(x, na.rm = TRUE)
)

res_agg <- res[, c(
  list(n_cpg = .N),
  lapply(.SD, agg_fun)
), by = .(region_id, chr, start, end), .SDcols = samples]

# filter
res_agg <- res_agg[n_cpg >= MIN_CPG]

fwrite(res_agg, OUTFILE, sep = "\t")
cat("Regions kept:", nrow(res_agg), "\n")
cat("Output:", OUTFILE, "\n")