
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