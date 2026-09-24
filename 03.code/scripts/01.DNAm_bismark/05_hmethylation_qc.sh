#!/bin/bash
# name£º05_hmethylation_qc.sh
# input£º*_bismark_bt2_pe.deduplicated.bismark.cov.gz
# output£º05_methylation_qc_hmeth.txt
# usage£º05_hmethylation_qc.sh <map_dir> <hmeth_dir> <outdir> <logdir>
# author:zhangyuxin,20260924

set -euo pipefail

MAP_DIR=$1
HMETH_DIR=$2
OUTDIR=$3
LOGDIR=$4

mkdir -p "$OUTDIR" "$LOGDIR"

LOG="${LOGDIR}/05_hmethylation_qc.log"
exec > >(tee -a "$LOG") 2>&1

HMAPP="${OUTDIR}/05_hmethylation_qc.txt"
PHAGE="NC_001416.1"

echo "Start hmethylation QC: $(date)"
echo "MAP_DIR:   $MAP_DIR"
echo "HMETH_DIR: $HMETH_DIR"
echo "OUTFILE:   $HMAPP"

echo -e "Sample\tConversion_efficiency(%)\tCpG_methylation(%)\tAvg_coverage" > "$HMAPP"

for dir in "${MAP_DIR}"/*/; do
    sample=$(basename "$dir")

    [[ "$sample" == *H* ]] || continue

    cov_gz=$(ls "${HMETH_DIR}/${sample}/"*deduplicated.bismark.cov.gz 2>/dev/null | head -1)

    if [[ -z "$cov_gz" ]]; then
        echo "WARNING: $sample cov.gz not found" >&2
        continue
    fi

    stats=$(zcat "$cov_gz" | awk -v phage="$PHAGE" '
        {
            met_all  += $5
            unmet_all += $6
            depth_all += ($5 + $6)
            n_all++
        }
        $1 == phage {
            met_phage   += $5
            unmet_phage += $6
        }
        END {
            total_phage = met_phage + unmet_phage
            conv_eff = (total_phage > 0) ? (unmet_phage/total_phage)*100 : "NA"

            total_all = met_all + unmet_all
            meth_all  = total_all > 0 ? (met_all/total_all)*100 : 0
            avg_cov   = n_all > 0 ? depth_all/n_all : "NA"

            printf "%.2f\t%.2f\t%.2f", conv_eff, meth_all, avg_cov
        }')

    echo -e "${sample}\t${stats}" >> "$HMAPP"
    echo "  H sample done: $sample"
done

echo "Hmethylation QC done: $(date)"
echo "  Rows: $(wc -l < "$HMAPP")"