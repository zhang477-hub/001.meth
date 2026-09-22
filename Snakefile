import pandas as pd

configfile: "03.code/config/config.yaml"

SAMPLES = pd.read_csv("03.code/config/samples.tsv", sep="\t")["sample_id"].tolist()

include: "03.code/rules/01_trim.smk"

rule all:
    input:
        expand(f"{config[paths][trimmed_dir]}/{{sample}}/{{sample}}_R1_val_1.fq.gz", sample = SAMPLES)