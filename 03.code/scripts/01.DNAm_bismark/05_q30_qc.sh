#!/bin/bash
# name£º05_q30.sh
# input£ºtrim before/after fastq.gz
# output£ºq30_stats.tsv
# usage£º05_q30.sh <raw_dir> <trim_dir> <outfile> <logdir> [r1_suffix] [r2_suffix] [tr1_suffix] [tr2_suffix]
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

LOG="${LOGDIR}/05_q30.log"
exec > >(tee -a "$LOG") 2>&1

echo "Start Q30 QC: $(date)"
echo "RAW_DIR:  $RAW_DIR"
echo "TRIM_DIR: $TRIM_DIR"
echo "OUTFILE:  $OUTFILE"
printf "sample\tq30_bases_before_b\tq30_rate_before(%%)\tq30_bases_after_b\tq30_rate_after(%%)\n" > "$OUTFILE"

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

    get_q30() {
        $SEQKIT stats -a -T -j "$THREADS" "$1" | awk -F'\t' 'NR==2 {
            # $5=sum_len, $15=Q30(%)
            printf "%d %.2f", $5, $15
        }'
    }

    read sum1_b q30r1_b <<< $(get_q30 "$before_r1")
    read sum2_b q30r2_b <<< $(get_q30 "$before_r2")
    read sum1_a q30r1_a <<< $(get_q30 "$after_r1")
    read sum2_a q30r2_a <<< $(get_q30 "$after_r2")

    total_bases_b=$((sum1_b + sum2_b))
    q30_bases_b=$(awk "BEGIN {printf \"%.0f\", $sum1_b * $q30r1_b / 100 + $sum2_b * $q30r2_b / 100}")
    q30_rate_b=$(awk "BEGIN {printf \"%.2f\", $total_bases_b ? $q30_bases_b / $total_bases_b * 100 : 0}")

    total_bases_a=$((sum1_a + sum2_a))
    q30_bases_a=$(awk "BEGIN {printf \"%.0f\", $sum1_a * $q30r1_a / 100 + $sum2_a * $q30r2_a / 100}")
    q30_rate_a=$(awk "BEGIN {printf \"%.2f\", $total_bases_a ? $q30_bases_a / $total_bases_a * 100 : 0}")

    printf "%s\t%s\t%s\t%s\t%s\n" \
        "$sample" "$q30_bases_b" "$q30_rate_b" "$q30_bases_a" "$q30_rate_a" >> "$OUTFILE"

    echo "  $sample done: Q30 before=${q30_rate_b}% after=${q30_rate_a}%"
done

echo "Q30 QC done: $(date)"
echo "  Rows: $(wc -l < "$OUTFILE")"