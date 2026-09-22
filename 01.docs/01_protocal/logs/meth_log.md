## 20260922

#### To do
- [ ] 写代码03.code/scripts/01.DNAm_bismark

#### code
- 03.code/scripts/01_qc/01_read_idat.R
- 03.code/rules/01_qc.smk

#### patra
- detection_pval = 0.01
- logFC_cutoff = 0.1

#### input
- 02.data/idat/
- 02.data/ref/sample_sheet.csv

#### output
- 04.results/qc/rgset.RData
- 04.results/qc/rgset_qc.pdf

#### figure/tables
- 04.results/figures/volcano.pdf ← 03.code/scripts/05_plot/volcano.R

#### tissue
- read.metharray.exp 报错，路径问题，已解决

#### Git
- commit: feat: 完成IDAT读取脚本

#### tommorow
- 写 normalize 脚本