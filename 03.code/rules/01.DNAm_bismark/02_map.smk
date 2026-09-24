rule bismark_map:
    input:
        r1 = f"{config['paths']['trim_dir']}/{{sample}}/{{sample}}_R1_val_1.fq.gz",
        r2 = f"{config['paths']['trim_dir']}/{{sample}}/{{sample}}_R2_val_2.fq.gz"
    output:
        bam = f"{config['paths']['map_dir']}/{{sample}}/{{sample}}_R1_val_1_bismark_bt2_pe.bam"
    params:
        genome = config["paths"]["genome"],
        logdir = f"{config['paths']['logs_dir']}/02_map"
    threads: config["map"]["threads"]
    resources:
        mem_mb = 64000,
        runtime = 2880
    conda:
        "/data/home/quj_lab/zhangyuxin/01-project/001.meth/03.code/envs/wgbs.yaml"
    shell:
        "source /data/home/quj_lab/zhangyuxin/anaconda3/etc/profile.d/conda.sh && conda run -n wgbs bash /data/home/quj_lab/zhangyuxin/01-project/001.meth/03.code/scripts/01.DNAm_bismark/02_map.sh {config[paths][trim_dir]} {wildcards.sample} {config[paths][map_dir]} {params.genome} {threads} {params.logdir}"