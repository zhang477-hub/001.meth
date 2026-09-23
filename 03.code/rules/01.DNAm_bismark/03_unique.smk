rule dedup_unique:
    input:
        bam = f"{config['paths']['map_dir']}/{{sample}}/{{sample}}_R1_val_1_bismark_bt2_pe.bam"
    output:
        bam = f"{config['paths']['unique_dir']}/{{sample}}/{{sample}}_R1_val_1_bismark_bt2_pe.deduplicated.bam",
        report = f"{config['paths']['unique_dir']}/{{sample}}/{{sample}}_R1_val_1_bismark_bt2_pe.deduplication_report.txt"
    params:
        logdir = f"{config['paths']['logs_dir']}/03_unique"
    threads: config["unique"]["threads"]
    resources:
        mem_mb = 64000,
        runtime = 2880
    conda:
        "03.code/envs/wgbs.yaml"
    shell:
        "bash /data/home/quj_lab/zhangyuxin/01-project/001.meth/03.code/scripts/01.DNAm_bismark/03_unique.sh "
        "{config[paths][map_dir]} {wildcards.sample} "
        "{config[paths][unique_dir]} "
        "{threads} {params.logdir}"
