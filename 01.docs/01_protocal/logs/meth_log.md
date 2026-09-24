## 20260923

#### To do
- [x] 写代码03.code/scripts/01.DNAm_bismark，整理01_trim.sh，02_map.sh，03_unique，04_meth.sh
- [x] 测试这些代码，输出在04.results/debug/bat/

#### code
- 03.code/scripts/01.DNAm_bismark/01_trim.sh
- 03.code/scripts/01.DNAm_bismark/02_map.sh
- 03.code/scripts/01.DNAm_bismark/03_unique.sh
- 03.code/scripts/01.DNAm_bismark/04_meth.sh
- 03.code/rules/01.DNAm_bismark/01_trim.smk
- 03.code/rules/01.DNAm_bismark/02_map.smk
- 03.code/rules/01.DNAm_bismark/03_unique.smk
- 03.code/rules/01.DNAm_bismark/04_meth.smk
- 03.code/config/01.DNAm_bismark_config.yaml

#### para
- trim: length: 50 quality: 30
- map：threads: 20
- unique： threads: 24
- meth：threads: 20，buffer_size: "10G"

#### Git
- commit: feat: 完成rules,config写作

## 20260924

#### To do
- [x] 测试 Snakemake dry run，写配置脚本
- [x] 写质控分析代码

#### code
- 03.code/scripts/run_snakemake_pbs.sh
- `03.code/scripts/01.DNAm_bismark/05_mapping_qc.sh`
- `03.code/scripts/01.DNAm_bismark/05_methylation_qc.sh`
- `03.code/scripts/01.DNAm_bismark/05_hmethylation_qc.sh`
- `03.code/rules/01.DNAm_bismark/05_qc.smk`



#### Git
- commit: feat: 完成QC脚本

