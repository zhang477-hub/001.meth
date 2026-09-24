#!/bin/bash
# usage£ºrun_snakemake_pbs.sh [¶îÍâ²ÎÊý]
# author£ºzhangyuxin, 20260924

set -euo pipefail

cd /data/home/quj_lab/zhangyuxin/01-project/001.meth

source /data/home/quj_lab/zhangyuxin/anaconda3/etc/profile.d/conda.sh

conda activate snakemake

snakemake \
  --cluster "qsub -l nodes=1:ppn={threads} -l mem={resources.mem_mb}mb -l walltime={resources.runtime}:00:00 -q batch -j oe -o /data/home/quj_lab/zhangyuxin/01-project/001.meth/04.results/01_DNAm_bismark/logs/pbs/" \
  --jobs 40 \
  --latency-wait 60 \
  --rerun-incomplete \
  --keep-going \
  "$@"