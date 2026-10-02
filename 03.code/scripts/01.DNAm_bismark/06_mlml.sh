#!/bin/bash
# name£º06_mlml.sh
# input£ºM sample cov.gz + matched H sample cov.gz
# output£º*_MLML_results.txt, *_check.txt
# usage£º06_mlml.sh <sample_m> <mc_dir> <hmc_dir> <outbase> <logdir>
# author:zhangyuxin,20260929

set -euo pipefail

SAMPLE_M=$1
MC_DIR=$2
HMC_DIR=$3
OUTBASE=$4
LOGDIR=${5:-${OUTBASE}/logs}

SAMPLE_H=${SAMPLE_M/M/H}

OUTDIR="${OUTBASE}/${SAMPLE_M}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

mkdir -p "$OUTDIR" "$LOGDIR"

LOG="${LOGDIR}/${SAMPLE_M}_mlml.log"
exec > >(tee -a "$LOG") 2>&1

echo "=================================================="
echo "Start 07_mlml: $(date)"
echo "SAMPLE_M: $SAMPLE_M"
echo "SAMPLE_H: $SAMPLE_H"
echo "MC_DIR:   $MC_DIR"
echo "HMC_DIR:  $HMC_DIR"
echo "OUTDIR:   $OUTDIR"
echo "=================================================="

mc_gz=$(ls "${MC_DIR}/${SAMPLE_M}/"*deduplicated.bismark.cov.gz 2>/dev/null | head -1)
hmc_gz=$(ls "${HMC_DIR}/${SAMPLE_H}/"*deduplicated.bismark.cov.gz 2>/dev/null | head -1)

if [[ -z "$mc_gz" ]]; then
    echo "ERROR: M cov.gz not found for $SAMPLE_M in $MC_DIR/$SAMPLE_M" >&2
    exit 1
fi
if [[ -z "$hmc_gz" ]]; then
    echo "ERROR: H cov.gz not found for $SAMPLE_H in $HMC_DIR/$SAMPLE_H" >&2
    exit 1
fi

echo "MC:  $mc_gz"
echo "HMC: $hmc_gz"

echo ""
echo "=== Step 1: 06_mlml.R ==="
Rscript "${SCRIPT_DIR}/06_mlml.R" \
  "$mc_gz" "$hmc_gz" "$SAMPLE_M" "$OUTDIR" "$LOGDIR"

MLML_OUT="${OUTDIR}/${SAMPLE_M}_MLML_results.txt"
[[ -f "$MLML_OUT" ]] || { echo "ERROR: $MLML_OUT not found" >&2; exit 1; }
echo "MLML results: $MLML_OUT"

echo ""
echo "=== Step 2: 06_check_mlml.R ==="
Rscript "${SCRIPT_DIR}/06_check_mlml.R" \
  "$MLML_OUT" "$SAMPLE_M" "$OUTDIR" "$LOGDIR" \
  "$mc_gz" "$hmc_gz"

CHECK_OUT="${OUTDIR}/${SAMPLE_M}_check.txt"
[[ -f "$CHECK_OUT" ]] || { echo "ERROR: $CHECK_OUT not found" >&2; exit 1; }
echo "Check results: $CHECK_OUT"

echo ""
echo "Sample $SAMPLE_M done: $(date)"