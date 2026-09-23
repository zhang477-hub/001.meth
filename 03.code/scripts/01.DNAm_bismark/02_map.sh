#!/bin/bash
# name£º02_map.sh
# input£º*_R1_val_1.fq.gz;*_R2_val_2.fq.gz
# output£º*_bismark_bt2_pe.bam
# usage£º02_map.sh <trim_dir> <sample> <outdir> <genome> <threads> <logdir>
# author:zhangyuxin,20260923

set -euo pipefail

TRIM_DIR=$1
SAMPLE=$2
OUTBASE=$3
GENOME=$4
THREADS=${5:-20}
LOGDIR=${6:-${OUTDIR}/logs}

OUTDIR="${OUTBASE}/${SAMPLE}"
mkdir -p "$OUTDIR" "$LOGDIR"

LOG="${LOGDIR}/${SAMPLE}_mapping.log"
exec > >(tee -a "$LOG") 2>&1

fq1=$(ls "${TRIM_DIR}/${SAMPLE}"/*_R1_val_1.fq.gz 2>/dev/null | head -1)
fq2=$(ls "${TRIM_DIR}/${SAMPLE}"/*_R2_val_2.fq.gz 2>/dev/null | head -1)

if [[ -z "$fq1" || -z "$fq2" ]]; then
    echo "ERROR: trimmed R1/R2 not found for $SAMPLE in $TRIM_DIR/$SAMPLE" >&2
    exit 1
fi

echo "Processing sample: $SAMPLE"
echo "R1: $fq1"
echo "R2: $fq2"
echo "Genome: $GENOME"
echo "Output: $OUTDIR"

bismark -p "$THREADS" -o "$OUTDIR" --non_directional --genome "$GENOME" -1 "$fq1" -2 "$fq2"

echo "Sample $SAMPLE done"