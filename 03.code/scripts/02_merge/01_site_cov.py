#!/usr/bin/env python3
# name：01_site_cov.py
# input：per sample, *_bismark.cov.gz
# output：beta_matrix.tsv, coverage_matrix.tsv, meth_counts.tsv, unmeth_counts.tsv
# usage：01_site_cov.py <meth_dir> <samples_file> <outdir> [--min-cov 5] [--mode intersect] [--chromosome chr1]
# author:zhangyuxin,20260929

import os
import sys
import gzip
import argparse

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("meth_dir")
    p.add_argument("samples_file")
    p.add_argument("outdir")
    p.add_argument("--min-cov", type=int, default=5)
    p.add_argument("--mode", choices=["intersect", "union"], default="intersect")
    p.add_argument("--chromosome", default=None)
    return p.parse_args()

def find_cov(meth_dir, sample):
    d = os.path.join(meth_dir, sample)
    if not os.path.isdir(d):
        return None
    for f in os.listdir(d):
        if f.endswith(".bismark.cov.gz"):
            return os.path.join(d, f)
    return None

def read_cov(path, min_cov, chrom=None):
    data = {}
    with gzip.open(path, "rt") as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 6:
                continue
            c, pos = parts[0], parts[1]
            if chrom and c != chrom:
                continue
            try:
                meth = int(parts[4]); unmeth = int(parts[5])
            except ValueError:
                continue
            if meth + unmeth < min_cov:
                continue
            data[(c, pos)] = (meth, unmeth)
    return data

def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    with open(args.samples_file) as f:
        samples = [l.strip() for l in f if l.strip() and not l.startswith("sample_id")]

    print(f"Samples: {len(samples)}", file=sys.stderr)

    sample_data = {}
    for i, s in enumerate(samples, 1):
        cov = find_cov(args.meth_dir, s)
        if cov is None:
            print(f"WARNING: no cov.gz for {s}", file=sys.stderr)
            continue
        print(f"[{i}/{len(samples)}] {s} ...", file=sys.stderr)
        sample_data[s] = read_cov(cov, args.min_cov, args.chromosome)
        print(f"  sites: {len(sample_data[s])}", file=sys.stderr)

    keysets = [set(d.keys()) for d in sample_data.values()]
    if args.mode == "intersect":
        common = set.intersection(*keysets)
    else:
        common = set.union(*keysets)

    print(f"Common sites: {len(common)}", file=sys.stderr)
    sorted_keys = sorted(common, key=lambda x: (x[0], int(x[1])))

    f_beta   = open(os.path.join(args.outdir, "beta_matrix.tsv"), "w")
    f_cov    = open(os.path.join(args.outdir, "coverage_matrix.tsv"), "w")
    f_meth   = open(os.path.join(args.outdir, "meth_counts.tsv"), "w")
    f_unmeth = open(os.path.join(args.outdir, "unmeth_counts.tsv"), "w")

    header = "chr\tpos\t" + "\t".join(samples) + "\n"
    for f in (f_beta, f_cov, f_meth, f_unmeth):
        f.write(header)

    for key in sorted_keys:
        c, p = key
        rows = {f_beta: [c, p], f_cov: [c, p], f_meth: [c, p], f_unmeth: [c, p]}
        for s in samples:
            if s not in sample_data or key not in sample_data[s]:
                rows[f_beta].append("NA"); rows[f_cov].append("0")
                rows[f_meth].append("0");  rows[f_unmeth].append("0")
                continue
            m, u = sample_data[s][key]
            total = m + u
            beta = m / total if total > 0 else 0
            rows[f_beta].append(f"{beta:.4f}")
            rows[f_cov].append(str(total))
            rows[f_meth].append(str(m))
            rows[f_unmeth].append(str(u))

        for f, row in rows.items():
            f.write("\t".join(row) + "\n")

    for f in (f_beta, f_cov, f_meth, f_unmeth):
        f.close()

    print("Done", file=sys.stderr)

if __name__ == "__main__":
    main()