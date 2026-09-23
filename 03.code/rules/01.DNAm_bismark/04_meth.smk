rule extract_meth:
    input:
        bam = f"{config['paths']['unique_dir']}/{{sample}}/{{sample}}_R1_val_1_bismark_bt2_pe.deduplicated.bam"
    output:
        cpg = f"{config['paths']['meth_dir']}/{{sample}}/{{sample}}_R1_val_1_bismark_bt2_pe.deduplicated.CpG_report.txt.gz"
    params:
        genome      = config["paths"]["genome"],
        buffer_size = config["meth"]["buffer_size"],
        logdir      = f"{config['paths']['logs_dir']}/04_meth"
    threads: config["meth"]["threads"]
    resources:
        mem_mb = 64000,
        runtime = 2880
    conda:
        "03.code/envs/wgbs.yaml"
    shell:
        "bash /data/home/quj_lab/zhangyuxin/01-project/001.meth/03.code/scripts/01.DNAm_bismark/04_meth.sh "
        "{config[paths][unique_dir]} {wildcards.sample} "
        "{config[paths][meth_dir]} "
        "{params.genome} {threads} {params.logdir}"