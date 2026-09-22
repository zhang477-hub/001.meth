#!/bin/bash
# name£º01_trim.sh
# input£º*R1_fq.gz;*R2_fq.gz
# output:*_R1_val_1.fq.gz;*_R2_val_2.fq.gz
# usage£º01_trim.sh <raw_dir> <sample> <outdir> <threads> <length> <quality> <logdir>
# author:zhangyuxin,20260922

set -euo pipefail

RAW_DIR=$1
SAMPLE=$2
OUTDIR=$3
THREADS=${4:-8}
LENGTH=${5:-50}
QUALITY=${6:-30}
LOGDIR=${7:-${OUTDIR}/logs}

mkdir -p "$OUTDIR" "$LOGDIR"

LOG="${LOGDIR}/${SAMPLE}_trim_galore.log"
exec > >(tee -a "$LOG") 2>&1

fq1=$(ls "${RAW_DIR}/${SAMPLE}"/*_R1.fq.gz 2>/dev/null | head -1)
fq2=$(ls "${RAW_DIR}/${SAMPLE}"/*_R2.fq.gz 2>/dev/null | head -1)

if [[ -z "$fq1" || -z "$fq2" ]]; then
    echo "ERROR: R1/R2 not found for $SAMPLE in $RAW_DIR/$SAMPLE" >&2
    exit 1
fi

echo "Processing sample: $SAMPLE"
echo "R1: $fq1"
echo "R2: $fq2"

trim_galore --paired --length "$LENGTH" -q "$QUALITY" --retain_unpaired --fastqc --gzip --cores "$THREADS" -o "$OUTDIR" "$fq1" "$fq2"

echo "Sample $SAMPLE done"