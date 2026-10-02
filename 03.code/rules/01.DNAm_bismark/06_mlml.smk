
rule run_mlml:
    input:
        mc_dir  = config["paths"]["meth_dir"],
        hmc_dir = config["paths"]["hmeth_dir"]
    output:
        mlml  = f"{config['paths']['mlml_dir']}/{{sample}}/{{sample}}_MLML_results.txt",
        check = f"{config['paths']['mlml_dir']}/{{sample}}/{{sample}}_check.txt"
    params:
        logdir = f"{config['paths']['logs_dir']}/06_mlml"
    threads: 1
    resources:
        mem_mb = 16000,
        runtime = 720
    shell:
        "source /data/home/quj_lab/zhangyuxin/anaconda3/etc/profile.d/conda.sh && "
        "conda run -n R4.4 bash "
        "/data/home/quj_lab/zhangyuxin/01-project/001.meth/03.code/scripts/01.DNAm_bismark/06_mlml.sh "
        "{wildcards.sample} "
        "{input.mc_dir} {input.hmc_dir} "
        "{config[paths][mlml_dir]} {params.logdir}"