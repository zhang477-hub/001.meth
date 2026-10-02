#!/bin/bash
# name£º04_sample_cov.sh
# input£ºper sample, 07_mlml + 05_qc
# output£ºglobal_meth.tsv, qc_all.tsv, mlml_check_all.tsv
# usage£º04_sample_cov.sh <mlml_dir> <qc_dir> <outdir> <samples_file> <logdir>
# author:zhangyuxin,20260929

set -euo pipefail

MLML_DIR=$1
QC_DIR=$2
OUTDIR=$3
SAMPLES_FILE=$4
LOGDIR=$5

mkdir -p "$OUTDIR" "$LOGDIR"

LOG="${LOGDIR}/04_sample_level.log"
exec > >(tee -a "$LOG") 2>&1

echo "Start 04_sample_level: $(date)"

# global meth
GLOBAL="${OUTDIR}/global_meth.tsv"
printf "sample\tglobal_mC\tglobal_hmC\tglobal_C\tn_sites\n" > "$GLOBAL"

while read -r sample; do
    [[ -z "$sample" ]] && continue
    [[ "$sample" == "sample_id" ]] && continue

    res="${MLML_DIR}/${sample}/${sample}_MLML_results.txt"
    if [[ ! -f "$res" ]]; then
        echo "WARNING: $res not found" >&2
        continue
    fi

    awk -v s="$sample" '
        NR>1 { mC+=$3; hmC+=$4; C+=$5; n++ }
        END {
            if (n>0) printf "%s\t%.4f\t%.4f\t%.4f\t%d\n", s, mC/n, hmC/n, C/n, n
        }' "$res" >> "$GLOBAL"
done < "$SAMPLES_FILE"

# mlml_check
CHECK="${OUTDIR}/mlml_check_all.tsv"
CHECK_FILES=$(find "$MLML_DIR" -name "*_check.txt" 2>/dev/null | sort)

if [[ -n "$CHECK_FILES" ]]; then
    head -1 $(echo "$CHECK_FILES" | head -1) > "$CHECK"
    for f in $CHECK_FILES; do
        tail -n +2 "$f" >> "$CHECK"
    done
fi

# qc_summary
QC_ALL="${OUTDIR}/qc_all.tsv"
QC_FILES=$(find "$QC_DIR" -name "*_qc.tsv" -o -name "*_stats.tsv" 2>/dev/null | sort)

if [[ -n "$QC_FILES" ]]; then
    head -1 $(echo "$QC_FILES" | head -1) > "$QC_ALL"
    for f in $QC_FILES; do
        tail -n +2 "$f" >> "$QC_ALL"
    done
fi

echo "Outputs:"
echo "  $GLOBAL"
echo "  $CHECK"
echo "  $QC_ALL"
echo "Done: $(date)"