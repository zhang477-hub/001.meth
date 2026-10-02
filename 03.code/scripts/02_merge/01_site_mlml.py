#!/usr/bin/env python3
# name：01_site_mlml.py
# input：per sample, *_MLML_results.txt
# output：mC_matrix.tsv, hmC_matrix.tsv, C_matrix.tsv
# usage：01_site_mlml.py <mlml_dir> <samples_file> <outdir> [--mode intersect] [--chromosome chr1]
# author:zhangyuxin,20260930

import os
import sys
import argparse

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("mlml_dir")
    p.add_argument("samples_file")
    p.add_argument("outdir")
    p.add_argument("--mode", choices=["intersect", "union"], default="intersect")
    p.add_argument("--chromosome", default=None)
    return p.parse_args()

def find_mlml(mlml_dir, sample):
    d = os.path.join(mlml_dir, sample)
    if not os.path.isdir(d):
        return None
    for f in os.listdir(d):
        if f.endswith("_MLML_results.txt"):
            return os.path.join(d, f)
    return None

def read_mlml(path, chrom=None):
    data = {}
    with open(path) as f:
        f.readline()  # header
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 5:
                continue
            c, pos = parts[0], parts[1]
            if chrom and c != chrom:
                continue
            try:
                mC  = float(parts[2])
                hmC = float(parts[3])
                C   = float(parts[4])
            except ValueError:
                continue
            data[(c, pos)] = (mC, hmC, C)
    return data

def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    with open(args.samples_file) as f:
        samples = [l.strip() for l in f
                   if l.strip() and not l.startswith("sample_id")]

    print(f"Samples: {len(samples)}", file=sys.stderr)

    sample_data = {}
    for i, s in enumerate(samples, 1):
        fpath = find_mlml(args.mlml_dir, s)
        if fpath is None:
            print(f"WARNING: no MLML_results for {s}", file=sys.stderr)
            continue
        print(f"[{i}/{len(samples)}] {s} ...", file=sys.stderr)
        sample_data[s] = read_mlml(fpath, args.chromosome)
        print(f"  sites: {len(sample_data[s])}", file=sys.stderr)

    if not sample_data:
        sys.exit("No data read")

    keysets = [set(d.keys()) for d in sample_data.values()]
    if args.mode == "intersect":
        common = set.intersection(*keysets)
    else:
        common = set.union(*keysets)

    print(f"Common sites: {len(common)}", file=sys.stderr)
    sorted_keys = sorted(common, key=lambda x: (x[0], int(x[1])))

    header = "chr\tpos\t" + "\t".join(samples) + "\n"

    f_mC  = open(os.path.join(args.outdir, "mC_matrix.tsv"),  "w")
    f_hmC = open(os.path.join(args.outdir, "hmC_matrix.tsv"), "w")
    f_C   = open(os.path.join(args.outdir, "C_matrix.tsv"),   "w")
    for f in (f_mC, f_hmC, f_C):
        f.write(header)

    for key in sorted_keys:
        c, p = key
        row_mC  = [c, p]
        row_hmC = [c, p]
        row_C   = [c, p]

        for s in samples:
            if s not in sample_data or key not in sample_data[s]:
                row_mC.append("NA"); row_hmC.append("NA"); row_C.append("NA")
                continue
            mC, hmC, C = sample_data[s][key]
            row_mC.append(f"{mC:.4f}")
            row_hmC.append(f"{hmC:.4f}")
            row_C.append(f"{C:.4f}")

        f_mC.write("\t".join(row_mC) + "\n")
        f_hmC.write("\t".join(row_hmC) + "\n")
        f_C.write("\t".join(row_C) + "\n")

    for f in (f_mC, f_hmC, f_C):
        f.close()

    print("Done", file=sys.stderr)
    print(f"  {args.outdir}/mC_matrix.tsv",  file=sys.stderr)
    print(f"  {args.outdir}/hmC_matrix.tsv", file=sys.stderr)
    print(f"  {args.outdir}/C_matrix.tsv",   file=sys.stderr)

if __name__ == "__main__":
    main()