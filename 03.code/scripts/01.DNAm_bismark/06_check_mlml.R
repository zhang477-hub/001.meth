#!/usr/bin/env Rscript
# name£º06_check_mlml.R
# input£º*_MLML_results.txt, ¿ÉÑ¡ *_bismark.cov.gz
# output£º*_check.txt
# usage£ºcheck_mlml.R <results_file> <sample> <outdir> <logdir> [wgbs_cov] [aceseq_cov]
# author:zhangyuxin,20260928

suppressPackageStartupMessages({
  library(MLML2R)
})

# arg
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 4) {
  stop("Usage: check_mlml.R <results_file> <sample> <outdir> <logdir> [wgbs_cov] [aceseq_cov]")
}

RESULTS_FILE <- args[1]
SAMPLE       <- args[2]
OUTDIR       <- args[3]
LOGDIR       <- args[4]
WGBS_FILE    <- if (length(args) >= 5) args[5] else NULL
ACESEQ_FILE  <- if (length(args) >= 6) args[6] else NULL

dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
dir.create(LOGDIR, recursive = TRUE, showWarnings = FALSE)

LOG_FILE <- file.path(LOGDIR, paste0(SAMPLE, "_check_mlml.log"))
OUTFILE  <- file.path(OUTDIR, paste0(SAMPLE, "_check.txt"))

log_con <- file(LOG_FILE, open = "at")
sink(log_con, append = TRUE)
sink(log_con, append = TRUE, type = "message")

cat("==================================================\n")
cat("Start check_mlml:", format(Sys.time()), "\n")
cat("Sample:       ", SAMPLE, "\n")
cat("Results file: ", RESULTS_FILE, "\n")
cat("WGBS cov:     ", ifelse(is.null(WGBS_FILE), "NA", WGBS_FILE), "\n")
cat("ACE-seq cov:  ", ifelse(is.null(ACESEQ_FILE), "NA", ACESEQ_FILE), "\n")
cat("OUTFILE:      ", OUTFILE, "\n")
cat("==================================================\n")

if (!file.exists(RESULTS_FILE)) {
  stop(paste("Results file not found:", RESULTS_FILE))
}

res <- read.table(RESULTS_FILE, header = TRUE, stringsAsFactors = FALSE)

required_cols <- c("chr", "pos", "mC", "hmC", "C")
if (!all(required_cols %in% colnames(res))) {
  stop("Results file must contain columns: chr, pos, mC, hmC, C")
}

cat("Number of sites:", nrow(res), "\n")


range_ok <- all(res$mC  >= 0 & res$mC  <= 1, na.rm = TRUE) &
            all(res$hmC >= 0 & res$hmC <= 1, na.rm = TRUE) &
            all(res$C   >= 0 & res$C   <= 1, na.rm = TRUE)

sum_vals <- res$mC + res$hmC + res$C
sum_ok   <- all(abs(sum_vals - 1) < 0.01, na.rm = TRUE)


mean_mC  <- mean(res$mC,  na.rm = TRUE)
mean_hmC <- mean(res$hmC, na.rm = TRUE)
mean_C   <- mean(res$C,   na.rm = TRUE)


cor_total <- NA
cor_hmC   <- NA
n_common  <- NA

if (!is.null(WGBS_FILE) && !is.null(ACESEQ_FILE)) {
  if (file.exists(WGBS_FILE) && file.exists(ACESEQ_FILE)) {

    read_cov <- function(file_path) {
      if (grepl("\\.gz$", file_path)) {
        con <- gzfile(file_path, "rt")
        on.exit(close(con))
      } else {
        con <- file_path
      }
      read.table(con, header = FALSE, stringsAsFactors = FALSE,
                 col.names = c("chr", "start", "end", "meth_perc",
                               "meth_reads", "unmeth_reads"))
    }

    wgbs   <- read_cov(WGBS_FILE)
    aceseq <- read_cov(ACESEQ_FILE)

    merged <- merge(wgbs[, c("chr", "start", "meth_reads", "unmeth_reads")],
                    aceseq[, c("chr", "start", "meth_reads", "unmeth_reads")],
                    by = c("chr", "start"),
                    suffixes = c("_wgbs", "_aceseq"))

    res2 <- res
    colnames(res2)[colnames(res2) == "pos"] <- "start"
    check_df <- merge(merged, res2, by = c("chr", "start"))
    n_common <- nrow(check_df)

    if (n_common > 0) {
      check_df$wgbs_total    <- check_df$meth_reads_wgbs + check_df$unmeth_reads_wgbs
      check_df$aceseq_total  <- check_df$meth_reads_aceseq + check_df$unmeth_reads_aceseq
      check_df$wgbs_ratio    <- check_df$meth_reads_wgbs / check_df$wgbs_total
      check_df$aceseq_ratio  <- check_df$meth_reads_aceseq / check_df$aceseq_total
      check_df$est_total     <- check_df$mC + check_df$hmC

      cor_total <- cor(check_df$est_total,  check_df$wgbs_ratio,   use = "complete.obs")
      cor_hmC   <- cor(check_df$hmC,        check_df$aceseq_ratio, use = "complete.obs")
    }
  } else {
    cat("WARNING: input cov file not found, skip consistency check\n")
  }
}


summary_df <- data.frame(
  sample           = SAMPLE,
  n_sites          = nrow(res),
  range_ok         = range_ok,
  sum_ok           = sum_ok,
  mean_mC          = round(mean_mC,  4),
  mean_hmC         = round(mean_hmC, 4),
  mean_C           = round(mean_C,   4),
  n_common         = ifelse(is.na(n_common),  NA, n_common),
  cor_total        = ifelse(is.na(cor_total), NA, round(cor_total, 4)),
  cor_hmC          = ifelse(is.na(cor_hmC),   NA, round(cor_hmC,   4))
)

write.table(summary_df, file = OUTFILE,
            sep = "\t", row.names = FALSE, quote = FALSE)


cat("Basic checks:\n")
cat("  range [0,1]:  ", ifelse(range_ok, "YES", "NO"), "\n")
cat("  sum ~ 1:      ", ifelse(sum_ok,   "YES", "NO"), "\n")
cat("  mean mC:      ", round(mean_mC,  4), "\n")
cat("  mean hmC:     ", round(mean_hmC, 4), "\n")
cat("  mean C:       ", round(mean_C,   4), "\n")

if (!is.na(cor_total)) {
  cat("Consistency (with input):\n")
  cat("  common sites: ", n_common, "\n")
  cat("  cor(mC+hmC ~ WGBS):   ", round(cor_total, 4), "\n")
  cat("  cor(hmC ~ ACE-seq):   ", round(cor_hmC,   4), "\n")
}

cat("Summary saved to:", OUTFILE, "\n")
cat("Sample", SAMPLE, "done\n")
cat("End check_mlml:", format(Sys.time()), "\n")

sink(type = "message")
sink()
close(log_con)