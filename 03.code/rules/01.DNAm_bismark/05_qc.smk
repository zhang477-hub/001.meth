
rule mapping_qc:
    input:
        map_dir = config["paths"]["map_dir"]
    output:
        mapp = f"{config['paths']['qc_dir']}/mapp.txt"
    params:
        logdir = f"{config['paths']['logs_dir']}/05_qc"
    threads: 1
    resources:
        mem_mb = 4000,
        runtime = 60
    shell:
        "bash /data/home/quj_lab/zhangyuxin/01-project/001.meth/03.code/scripts/01.DNAm_bismark/05_mapping_qc.sh "
        "{input.map_dir} {output.mapp} {params.logdir}"




rule methylation_qc:
    input:
        map_dir  = config["paths"]["map_dir"],
        meth_dir = config["paths"]["meth_dir"]
    output:
        mapp = f"{config['paths']['qc_dir']}/05_methylation_qc.txt"
    params:
        logdir = f"{config['paths']['logs_dir']}/05_qc"
    threads: 1
    resources:
        mem_mb = 8000,
        runtime = 120
    shell:
        "bash /data/home/quj_lab/zhangyuxin/01-project/001.meth/03.code/scripts/01.DNAm_bismark/05_methylation_qc.sh "
        "{input.map_dir} {input.meth_dir} "
        "{config[paths][qc_dir]} {params.logdir}"




rule hmethylation_qc:
    input:
        map_dir   = config["paths"]["map_dir"],
        hmeth_dir = config["paths"]["hmeth_dir"]
    output:
        hmapp = f"{config['paths']['qc_dir']}/05_hmethylation_qc.txt"
    params:
        logdir = f"{config['paths']['logs_dir']}/05_qc"
    threads: 1
    resources:
        mem_mb = 8000,
        runtime = 120
    shell:
        "bash /data/home/quj_lab/zhangyuxin/01-project/001.meth/03.code/scripts/01.DNAm_bismark/05_hmethylation_qc.sh "
        "{input.map_dir} {input.hmeth_dir} "
        "{config[paths][qc_dir]} {params.logdir}"
        

rule q30_qc:
    input:
        raw_dir  = config["paths"]["raw_dir"],
        trim_dir = config["paths"]["trim_dir"]
    output:
        tsv = f"{config['paths']['qc_dir']}/q30_stats.tsv"
    params:
        logdir     = f"{config['paths']['logs_dir']}/05_qc",
        r1_suffix  = config["qc"]["r1_suffix"],
        r2_suffix  = config["qc"]["r2_suffix"],
        tr1_suffix = config["qc"]["trim_r1_suffix"],
        tr2_suffix = config["qc"]["trim_r2_suffix"]
    threads: 4
    resources:
        mem_mb = 4000,
        runtime = 120
    shell:
        "source /data/home/quj_lab/zhangyuxin/anaconda3/etc/profile.d/conda.sh && "
        "conda run -n qc bash /data/home/quj_lab/zhangyuxin/01-project/001.meth/03.code/scripts/01.DNAm_bismark/05_q30.sh "
        "{input.raw_dir} {input.trim_dir} {output.tsv} {params.logdir} "
        "{params.r1_suffix} {params.r2_suffix} "
        "{params.tr1_suffix} {params.tr2_suffix}"




rule fastq_qc:
    input:
        raw_dir  = config["paths"]["raw_dir"],
        trim_dir = config["paths"]["trim_dir"]
    output:
        tsv = f"{config['paths']['qc_dir']}/qc_stats.tsv"
    params:
        logdir     = f"{config['paths']['logs_dir']}/05_qc",
        r1_suffix  = config["qc"]["r1_suffix"],
        r2_suffix  = config["qc"]["r2_suffix"],
        tr1_suffix = config["qc"]["trim_r1_suffix"],
        tr2_suffix = config["qc"]["trim_r2_suffix"]
    threads: 4
    resources:
        mem_mb = 4000,
        runtime = 300
    shell:
        "source /data/home/quj_lab/zhangyuxin/anaconda3/etc/profile.d/conda.sh && "
        "conda run -n qc bash /data/home/quj_lab/zhangyuxin/01-project/001.meth/03.code/scripts/01.DNAm_bismark/05_fastq_qc.sh "
        "{input.raw_dir} {input.trim_dir} {output.tsv} {params.logdir} "
        "{params.r1_suffix} {params.r2_suffix} "
        "{params.tr1_suffix} {params.tr2_suffix}"        