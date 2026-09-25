#!/bin/bash
# name£º05_fastq_qc.sh
# input£ºtrim before/after fastq.gz
# output£ºqc_stats.tsv
# usage£º05_fastq_qc.sh <raw_dir> <trim_dir> <outfile> <logdir> [r1_suffix] [r2_suffix] [tr1_suffix] [tr2_suffix]
# author:zhangyuxin,20260925

set -euo pipefail

RAW_DIR=$1
TRIM_DIR=$2
OUTFILE=$3
LOGDIR=$4
R1_SUFFIX=${5:-_1.fq.gz}
R2_SUFFIX=${6:-_2.fq.gz}
TRIM_R1_SUFFIX=${7:-_1_val_1.fq.gz}
TRIM_R2_SUFFIX=${8:-_2_val_2.fq.gz}

SEQKIT="seqkit"
THREADS=4

mkdir -p "$(dirname "$OUTFILE")" "$LOGDIR"

LOG="${LOGDIR}/05_fastq_qc.log"
exec > >(tee -a "$LOG") 2>&1

echo "Start FASTQ QC: $(date)"
echo "RAW_DIR:  $RAW_DIR"
echo "TRIM_DIR: $TRIM_DIR"
echo "OUTFILE:  $OUTFILE"

printf "sample\ttotal_reads_b\ttotal_bases_b\tq20_bases_b\tq30_bases_b\tmean_length_b\tgc_b\ttotal_reads_a\ttotal_bases_a\tq20_bases_a\tq30_bases_a\tmean_length_a\tgc_a\n" > "$OUTFILE"

for sample_dir in "$TRIM_DIR"/*/; do
    sample=$(basename "$sample_dir")

    before_r1="$RAW_DIR/$sample/${sample}${R1_SUFFIX}"
    before_r2="$RAW_DIR/$sample/${sample}${R2_SUFFIX}"
    after_r1="$TRIM_DIR/$sample/${sample}${TRIM_R1_SUFFIX}"
    after_r2="$TRIM_DIR/$sample/${sample}${TRIM_R2_SUFFIX}"

    if [[ ! -f "$after_r1" || ! -f "$after_r2" ]]; then
        echo "WARNING: trimmed files missing for $sample" >&2
        continue
    fi
    if [[ ! -f "$before_r1" || ! -f "$before_r2" ]]; then
        echo "WARNING: raw files missing for $sample" >&2
        continue
    fi

    get_stats() {
        $SEQKIT stats -a -T -j "$THREADS" "$1" | awk -F'\t' 'NR==2 {
            printf "%d %d %.2f %.2f %.2f", $4, $5, $14, $15, $17
        }'
    }

    read num1_b sum1_b q20r1_b q30r1_b gc1_b <<< $(get_stats "$before_r1")
    read num2_b sum2_b q20r2_b q30r2_b gc2_b <<< $(get_stats "$before_r2")
    read num1_a sum1_a q20r1_a q30r1_a gc1_a <<< $(get_stats "$after_r1")
    read num2_a sum2_a q20r2_a q30r2_a gc2_a <<< $(get_stats "$after_r2")

    total_reads_b=$((num1_b + num2_b))
    total_bases_b=$((sum1_b + sum2_b))
    mean_length_b=$(awk "BEGIN {printf \"%.2f\", $total_reads_b ? $total_bases_b / $total_reads_b : 0}")
    q20_bases_b=$(awk "BEGIN {printf \"%.0f\", $sum1_b * $q20r1_b / 100 + $sum2_b * $q20r2_b / 100}")
    q30_bases_b=$(awk "BEGIN {printf \"%.0f\", $sum1_b * $q30r1_b / 100 + $sum2_b * $q30r2_b / 100}")
    gc_b=$(awk "BEGIN {printf \"%.2f\", $total_bases_b ? ($sum1_b * $gc1_b + $sum2_b * $gc2_b) / $total_bases_b : 0}")

    total_reads_a=$((num1_a + num2_a))
    total_bases_a=$((sum1_a + sum2_a))
    mean_length_a=$(awk "BEGIN {printf \"%.2f\", $total_reads_a ? $total_bases_a / $total_reads_a : 0}")
    q20_bases_a=$(awk "BEGIN {printf \"%.0f\", $sum1_a * $q20r1_a / 100 + $sum2_a * $q20r2_a / 100}")
    q30_bases_a=$(awk "BEGIN {printf \"%.0f\", $sum1_a * $q30r1_a / 100 + $sum2_a * $q30r2_a / 100}")
    gc_a=$(awk "BEGIN {printf \"%.2f\", $total_bases_a ? ($sum1_a * $gc1_a + $sum2_a * $gc2_a) / $total_bases_a : 0}")

    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
        "$sample" \
        "$total_reads_b" "$total_bases_b" "$q20_bases_b" "$q30_bases_b" "$mean_length_b" "$gc_b" \
        "$total_reads_a" "$total_bases_a" "$q20_bases_a" "$q30_bases_a" "$mean_length_a" "$gc_a" \
        >> "$OUTFILE"

    echo "  $sample done"
done

echo "FASTQ QC done: $(date)"
echo "  Rows: $(wc -l < "$OUTFILE")"