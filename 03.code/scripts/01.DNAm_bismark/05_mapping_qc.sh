#!/bin/bash
# name£º05_mapping_qc.sh
# input£º*_bismark_bt2_PE_report.txt
# output£ºmapp.txt
# usage£º05_mapping_qc.sh <map_dir> <outfile> <logdir>
# author:zhangyuxin,20260924

set -euo pipefail

MAP_DIR=$1
OUTFILE=$2
LOGDIR=$3

mkdir -p "$(dirname "$OUTFILE")" "$LOGDIR"

LOG="${LOGDIR}/05_mapping_qc.log"
exec > >(tee -a "$LOG") 2>&1

echo "Start mapping QC: $(date)"
echo "MAP_DIR: $MAP_DIR"
echo "OUTFILE: $OUTFILE"

{
    printf "sample\tTotal\tUnique\tMapping_efficiency\n"
    for dir in "${MAP_DIR}"/*/; do
        sample=$(basename "$dir")
        report=$(ls "${dir}"/*_bismark_bt2_PE_report.txt 2>/dev/null | head -1)

        if [[ -z "$report" ]]; then
            echo "WARNING: report not found for $sample" >&2
            continue
        fi

        awk -v sample="$sample" '
            /Sequence pairs analysed in total:/ {total=$NF}
            /Number of paired-end alignments with a unique best hit:/ {unique=$NF}
            /Mapping efficiency:/ {eff=$NF}
            END {
                if (total != "" && unique != "" && eff != "") {
                    printf "%s\t%s\t%s\t%s\n", sample, total, unique, eff
                }
            }
        ' "$report"
    done
} > "$OUTFILE"

echo "Mapping QC done: $(date)"
echo "Rows: $(wc -l < "$OUTFILE")"