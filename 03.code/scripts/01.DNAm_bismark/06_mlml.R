#!/usr/bin/env Rscript
# name£º06_mlml.R
# input£º*_bismark_bt2_pe.deduplicated.bismark.cov.gz (WGBS + ACE-seq)
# output£º*_MLML_results.txt
# usage£ºrun_mlml2.R <wgbs_cov> <aceseq_cov> <sample> <outdir> <logdir>
# author:zhangyuxin,20260928

suppressPackageStartupMessages({
  library(MLML2R)
})

# arg
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 5) {
  stop("Usage: run_mlml2.R <wgbs_cov> <aceseq_cov> <sample> <outdir> <logdir>")
}

WGBS_FILE   <- args[1]
ACESEQ_FILE <- args[2]
SAMPLE      <- args[3]
OUTDIR      <- args[4]
LOGDIR      <- args[5]


dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
dir.create(LOGDIR, recursive = TRUE, showWarnings = FALSE)

LOG_FILE <- file.path(LOGDIR, paste0(SAMPLE, "_mlml.log"))
OUTFILE  <- file.path(OUTDIR, paste0(SAMPLE, "_MLML_results.txt"))


log_con <- file(LOG_FILE, open = "at")
sink(log_con, append = TRUE)
sink(log_con, append = TRUE, type = "message")

cat("==================================================\n")
cat("Start run_mlml2:", format(Sys.time()), "\n")
cat("Sample:     ", SAMPLE, "\n")
cat("WGBS:       ", WGBS_FILE, "\n")
cat("ACE-seq:    ", ACESEQ_FILE, "\n")
cat("OUTFILE:    ", OUTFILE, "\n")
cat("==================================================\n")


if (!file.exists(WGBS_FILE)) {
  stop(paste("WGBS file not found:", WGBS_FILE))
}
if (!file.exists(ACESEQ_FILE)) {
  stop(paste("ACE-seq file not found:", ACESEQ_FILE))
}


read_cov <- function(file_path) {
  if (grepl("\\.gz$", file_path)) {
    con <- gzfile(file_path, "rt")
    on.exit(close(con))
  } else {
    con <- file_path
  }
  dat <- read.table(con, header = FALSE, stringsAsFactors = FALSE,
                    col.names = c("chr", "start", "end", "meth_perc",
                                  "meth_reads", "unmeth_reads"))
  return(dat)
}

cat("Reading WGBS file:", WGBS_FILE, "\n")
wgbs <- read_cov(WGBS_FILE)

cat("Reading ACE-seq file:", ACESEQ_FILE, "\n")
aceseq <- read_cov(ACESEQ_FILE)

cat("Merging datasets...\n")
merged <- merge(wgbs[, c("chr", "start", "meth_reads", "unmeth_reads")],
                aceseq[, c("chr", "start", "meth_reads", "unmeth_reads")],
                by = c("chr", "start"),
                suffixes = c("_wgbs", "_aceseq"))

if (nrow(merged) == 0) {
  stop("No common CpG sites after merging. Check input files.")
}

cat("Number of common CpG sites:", nrow(merged), "\n")

wgbs_matrix   <- as.matrix(merged[, c("meth_reads_wgbs",   "unmeth_reads_wgbs")])
aceseq_matrix <- as.matrix(merged[, c("meth_reads_aceseq", "unmeth_reads_aceseq")])

wgbs_meth   <- matrix(wgbs_matrix[, 1],   ncol = 1)   # C reads (5mC + 5hmC)
wgbs_unmeth <- matrix(wgbs_matrix[, 2],   ncol = 1)   # T reads (unmodified C)

aceseq_meth   <- matrix(aceseq_matrix[, 1], ncol = 1) # C reads (5hmC)
aceseq_unmeth <- matrix(aceseq_matrix[, 2], ncol = 1) # T reads (5mC + unmodified C)

cat("Running MLML estimation...\n")

result <- MLML(
  T.matrix  = wgbs_meth,
  U.matrix  = wgbs_unmeth,
  G.matrix  = aceseq_unmeth,
  H.matrix  = aceseq_meth,
  iterative = FALSE,
  tol       = 1e-5
)

output <- data.frame(
  chr = merged$chr,
  pos = merged$start,
  mC  = result$mC,
  hmC = result$hmC,
  C   = result$C
)

valid <- (output$mC  >= 0) & (output$hmC >= 0) & (output$C >= 0) &
         (output$mC + output$hmC + output$C <= 1 + 1e-6)
output_filtered <- output[valid, ]

cat("Sites before filtering:", nrow(output), "\n")
cat("Sites after filtering:", nrow(output_filtered), "\n")

write.table(output_filtered, file = OUTFILE,
            sep = "\t", row.names = FALSE, quote = FALSE)

cat("Results saved to:", OUTFILE, "\n")
cat("Sample", SAMPLE, "done\n")
cat("End run_mlml2:", format(Sys.time()), "\n")

sink(type = "message")
sink()
close(log_con)