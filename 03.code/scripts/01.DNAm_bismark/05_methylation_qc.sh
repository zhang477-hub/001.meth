#!/bin/bash
# name£º05_methylation_qc.sh
# input£º*_bismark_bt2_pe.deduplicated.bismark.cov.gz, Non_CpG_context_*.txt.gz
# output£º05_methylation_qc.txt
# usage£º05_methylation_qc.sh <map_dir> <meth_dir> <outdir> <logdir>
# author:zhangyuxin,20260924

set -euo pipefail

MAP_DIR=$1
METH_DIR=$2
OUTDIR=$3
LOGDIR=$4

mkdir -p "$OUTDIR" "$LOGDIR"

LOG="${LOGDIR}/05_methylation_qc.log"
exec > >(tee -a "$LOG") 2>&1

MAPP="${OUTDIR}/05_methylation_qc.txt"
PHAGE="NC_001416.1"

echo "Start methylation QC: $(date)"
echo "MAP_DIR:  $MAP_DIR"
echo "METH_DIR: $METH_DIR"
echo "OUTFILE:  $MAPP"

echo -e "Sample\tConversion_efficiency(%)\tCpG_methylation(%)\tAvg_coverage" > "$MAPP"

for dir in "${MAP_DIR}"/*/; do
    sample=$(basename "$dir")

    [[ "$sample" == *M* ]] || continue

    noncpg_gz=$(ls "${METH_DIR}/${sample}/Non_CpG_context_"*deduplicated.txt.gz 2>/dev/null | head -1)
    cov_gz=$(ls "${METH_DIR}/${sample}/"*deduplicated.bismark.cov.gz 2>/dev/null | head -1)

    if [[ -z "$noncpg_gz" || -z "$cov_gz" ]]; then
        echo "WARNING: missing files for $sample" >&2
        continue
    fi

    conv_eff=$(zcat "$noncpg_gz" | awk -v phage="$PHAGE" '
        NR==1 {next}
        $3 == phage {
            s=$NF
            if (s ~ /[XH]/) met++
            else if (s ~ /[xh]/) unmet++
        }
        END {
            total = met + unmet
            printf "%.2f", total ? (unmet/total)*100 : "NA"
        }')

    cov_stats=$(zcat "$cov_gz" | awk '
        { met+=$5; unmet+=$6; depth+=$5+$6; n++ }
        END {
            total = met + unmet
            meth  = total ? (met/total)*100 : 0
            avg   = n ? depth/n : 0
            printf "%.2f\t%.2f", meth, avg
        }')

    echo -e "${sample}\t${conv_eff}\t${cov_stats}" >> "$MAPP"
    echo "  M sample done: $sample"
done

echo "Methylation QC done: $(date)"
echo "  Rows: $(wc -l < "$MAPP")"