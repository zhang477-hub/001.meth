#!/bin/bash
# name£º03_unique.sh
# input£º*_bismark_bt2_pe.bam
# output£º*_bismark_bt2_pe.deduplicated.bam
# usage£º03_unique.sh <map_dir> <sample> <outdir> <threads> <logdir>
# author:zhangyuxin,20260923

set -euo pipefail

MAP_DIR=$1
SAMPLE=$2
OUTBASE=$3
THREADS=${4:-24}
LOGDIR=${5:-${OUTDIR}/logs}

OUTDIR="${OUTBASE}/${SAMPLE}"
mkdir -p "$OUTDIR" "$LOGDIR"

LOG="${LOGDIR}/${SAMPLE}_unique.log"
exec > >(tee -a "$LOG") 2>&1

bam=$(ls "${MAP_DIR}/${SAMPLE}"/*_bismark_bt2_pe.bam 2>/dev/null | head -1)

if [[ -z "$bam" ]]; then
    echo "ERROR: bam not found for $SAMPLE in $MAP_DIR/$SAMPLE" >&2
    exit 1
fi

echo "Processing sample: $SAMPLE"
echo "BAM: $bam"
echo "Output: $OUTDIR"
echo "Threads: $THREADS"

# 1. dedup
deduplicate_bismark --bam "$bam" --output_dir "$OUTDIR" --paired
dedup_bam=$(ls "${OUTDIR}"/*.deduplicated.bam 2>/dev/null | head -1)

if [[ -z "$dedup_bam" ]]; then
    echo "ERROR: deduplicated bam not found in $OUTDIR" >&2
    exit 1
fi

# 2. sort
#srt_bam="${OUTDIR}/${SAMPLE}.sort.rmdup.bam"
#samtools sort -@ "$THREADS" -l 9 -o "$srt_bam" "$dedup_bam"

# 3. index
#samtools index -@ "$THREADS" "$srt_bam"

# 4. statics
samtools flagstat -@ "$THREADS" "$dedup_bam" > "${OUTDIR}/${SAMPLE}.sort.rmdup.flagstat.txt"
#samtools idxstats -@ "$THREADS" "$srt_bam" > "${OUTDIR}/${SAMPLE}.sort.rmdup.idxstats.txt"

echo "Sample $SAMPLE done"