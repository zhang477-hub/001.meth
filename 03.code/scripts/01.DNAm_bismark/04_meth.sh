#!/bin/bash
# name£º04_meth.sh
# input£º*_bismark_bt2_pe.deduplicated.bam
# output£º*_CpG_report.txt.gz;*.bedGraph.gz
# usage£º04_meth.sh <unique_dir> <sample> <outdir> <genome> <threads> <logdir>
# author:zhangyuxin,20260923

set -euo pipefail

UNIQUE_DIR=$1
SAMPLE=$2
OUTBASE=$3
GENOME=$4
THREADS=${5:-20}
LOGDIR=${6:-${OUTDIR}/logs}

OUTDIR="${OUTBASE}/${SAMPLE}"
mkdir -p "$OUTDIR" "$LOGDIR"

LOG="${LOGDIR}/${SAMPLE}_meth.log"
exec > >(tee -a "$LOG") 2>&1

bam=$(ls "${UNIQUE_DIR}/${SAMPLE}"/*.deduplicated.bam 2>/dev/null | head -1)

if [[ -z "$bam" ]]; then
    echo "ERROR: deduplicated bam not found for $SAMPLE in $UNIQUE_DIR/$SAMPLE" >&2
    exit 1
fi

echo "Processing sample: $SAMPLE"
echo "BAM: $bam"
echo "Genome: $GENOME"
echo "Output: $OUTDIR"
echo "Threads: $THREADS"

bismark_methylation_extractor -p --gzip --bedGraph --buffer_size 10G --parallel "$THREADS" --cytosine_report --comprehensive --counts --merge_non_CpG --genome_folder "$GENOME" "$bam" -o "$OUTDIR"

echo "Sample $SAMPLE done"