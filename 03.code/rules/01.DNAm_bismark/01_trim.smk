rule trim_galore:
    input:
        raw_dir = lambda w: f"{config['paths']['raw_dir']}/{w.sample}"
    output:
        r1 = f"{config['paths']['trim_dir']}/{{sample}}/{{sample}}_R1_val_1.fq.gz",
        r2 = f"{config['paths']['trim_dir']}/{{sample}}/{{sample}}_R2_val_2.fq.gz"
    params:
        length  = config["trim"]["length"],
        quality = config["trim"]["quality"],
        logdir  = f"{config['paths']['logs_dir']}/01_trim"
    threads: 1
    resources:
        mem_mb = 8000,
        runtime = 1440
    shell:
        "source /data/home/quj_lab/zhangyuxin/anaconda3/etc/profile.d/conda.sh && conda run -n trim bash /data/home/quj_lab/zhangyuxin/01-project/001.meth/03.code/scripts/01.DNAm_bismark/01_trim.sh {input.raw_dir} {wildcards.sample} {config[paths][trim_dir]} {threads} {params.length} {params.quality} {params.logdir}"