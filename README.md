# RNA-seq 分析流水线（R / DESeq2）

一套分步运行的 RNA-seq 分析代码，覆盖从原始 count 矩阵到差异分析、靶标筛选、功能富集、可视化和自动报告的完整流程。

本仓库**只包含分析代码**，不包含表达矩阵、样本信息、注释表、富集结果、图片或研究结论。

## 流程结构

```
Pipeline/
├── 00_config.R          公共配置、依赖检查、路径与样本校验
├── 01_read_data.R       读取 count 矩阵与注释
├── 02_filter.R          低表达过滤与注释 1:1 对齐
├── 03_qc.R              vst + PCA 样本质控
├── 04_deseq.R           DESeq2 差异分析 + lfcShrink
├── 05_targets.R         过表达/敲低反向一致性靶标筛选
├── 06_go.R              GO 注释入口（本地文件或 UniProt 回退）
├── 07_enrichment.R      GO / KEGG / PFAM / TF 家族富集
├── 08_visualization.R   火山图、维恩图、热图、富集气泡图
├── 09_summary.R         关键统计与汇总表
├── 10_report.R          由中间结果自动生成 Markdown 报告
├── build_gene_go.R      参考基因组 GFF + eggNOG 注释桥接到 gene_id
└── run_eggnog.sh        离线运行 eggNOG-mapper
```

## 环境依赖

- R >= 4.3
- CRAN：`ggplot2`、`pheatmap`、`RColorBrewer`、`svglite`
- Bioconductor：`DESeq2`、`limma`、`apeglm`

`00_config.R` 会检查并尝试安装缺失的 R 包。也可以手动安装：

```r
install.packages(c("ggplot2", "pheatmap", "RColorBrewer", "svglite"))
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("DESeq2", "limma", "apeglm"))
```

## 输入数据

默认读取一个制表符分隔的 count 矩阵：

```
gene_id    sample_1    sample_2    sample_3    ...
gene_1     120         98          143         ...
gene_2     0           2           1           ...
```

要求：

- 第一列为 `gene_id`，且不能重复；
- 其余列为整数 read count；
- 下游按 `sample_cols` 中的**列名**取样本，不依赖文件列序。

把输入文件放到项目目录下，或通过环境变量指定：

```bash
export RNASEQ_BASE=/path/to/project
export RNASEQ_COUNT=/path/to/gene_count.xls
```

## 配置样本与分组

编辑 `Pipeline/00_config.R` 中的 `sample_cols` 和 `cond`：

```r
sample_cols <- c("ctrl_1", "ctrl_2",
                 "oe_1",   "oe_2",
                 "kd_1",   "kd_2")

cond <- factor(rep(c("ctrl", "oe", "kd"), each = 2),
               levels = c("ctrl", "oe", "kd"))
```

`cond` 的顺序必须与 `sample_cols` 一一对应。参考组建议放在 `levels` 的第一个位置。

> 仓库中的样本名和分组仅为示意。不同项目应替换为自己的样本列与实验设计。

## 运行

按顺序运行，任何一步失败都会停在上游：

```bash
for s in 01_read_data 02_filter 03_qc \
         04_deseq 05_targets 06_go \
         07_enrichment 08_visualization \
         09_summary 10_report; do
  Rscript Pipeline/$s.R || break
done
```

也可以单独重跑某一步。所有中间结果保存在 `results/rds/`，表格和图分别写入 `results/tables/` 与 `results/figures/`。

## 常用参数

| 环境变量 | 作用 | 默认值 |
|---|---|---|
| `RNASEQ_BASE` | 项目根目录 | `~/rna-seq` |
| `RNASEQ_COUNT` | count 矩阵路径 | `$RNASEQ_BASE/gene_count.xls` |
| `FILTER_MIN_COUNT` | 最低 read 计数 | `10` |
| `FILTER_MIN_SAMPLES` | 至少多少个样本达标 | `2` |
| `DEG_PADJ` | FDR 阈值 | `0.05` |
| `DEG_LFC` | `\|log2FC\|` 阈值 | `1` |
| `DESEQ_BATCH` | 与样本同序的批次向量 | 空 |
| `DESEQ_SHRINK` | `apeglm` 或 `normal` | `apeglm` |
| `PCA_NTOP` | PCA 使用的高变基因数 | `500` |
| `GO_ANNOT_FILE` | 本地 gene -> GO 注释文件 | 空 |
| `KEGG_DIR` | KO / pathway 定义文件目录 | `~/eggnog_data` |
| `TARGET_GENE_SYMBOL` | 报告中用于按名检索的目标基因符号 | `TF` |

示例：

```bash
FILTER_MIN_COUNT=10 \
FILTER_MIN_SAMPLES=2 \
DEG_PADJ=0.05 \
DEG_LFC=1 \
Rscript Pipeline/04_deseq.R
```

## 功能注释

### 方式一：本地 gene -> GO 文件

提供一个带表头的文件，至少包含 gene ID 与 GO ID 两列。`06_go.R` 会自动识别常见列名，并支持一个单元格内用逗号、分号或竖线分隔多个 GO：

```bash
GO_ANNOT_FILE=/path/to/go_by_gene.tsv \
Rscript Pipeline/06_go.R
```

### 方式二：参考基因组 GFF + eggNOG

1. 对蛋白 FASTA 运行 eggNOG-mapper：

```bash
PROTEIN_FASTA=/path/to/protein.faa \
EGGNOG_DATA=~/eggnog_data \
EGGNOG_CPU=8 \
bash Pipeline/run_eggnog.sh
```

2. 用 GFF 把蛋白 ID 映射回基因 ID，并生成 GO / KO 文件：

```bash
GFF_FILE=/path/to/genomic.gff \
PROT_GO_FILE=~/eggnog_data/emapper/out.emapper.annotations \
GO_OBO=~/eggnog_data/go-basic.obo \
Rscript Pipeline/build_gene_go.R
```

3. 将生成的文件交给 06 / 07：

```bash
GO_ANNOT_FILE=~/eggnog_data/go_by_gene.tsv \
Rscript Pipeline/06_go.R
Rscript Pipeline/07_enrichment.R
```

## 输出

```
results/
├── rds/          各步骤中间对象（可单独重跑）
├── tables/       DEG、靶标、富集、QC 坐标、方法说明
├── figures/      PCA、火山图、热图、维恩图、富集气泡图
├── summary.csv
└── ANALYSIS_REPORT.md
```

关键方法说明会随流程落盘，包括预处理、差异分析、靶标筛选、富集背景和可视化参数。

## 几个实现细节

- **样本按名对齐**：只校验列是否存在，下游按 `sample_cols` 重排，避免文件列序变化导致分组错位。
- **注释按 gene_id 对齐**：过滤后用 `match()` 重建 1:1 对应关系，出现 NA 立即停止。
- **先 QC 再 DE**：PCA 使用 `vst(blind = TRUE)`，不提前使用分组信息。
- **LFC 收缩口径**：`apeglm` 只接受 `coef`；两个处理组之间的对比使用两个后验均值之差，没有可用的后验 SE，不能用收缩后的 LFC 反推 z 值。
- **富集背景**：使用“有机会被判为 DEG”的基因（`padj` 非 NA），而不是全部过滤基因；旧口径保留为敏感性对照。
- **注释来源留痕**：记录本地注释文件、匹配率、GO 密度；外部 KEGG 文件记录路径、大小、md5 和修改时间。
- **报告现算**：`09_summary.R` 与 `10_report.R` 从 rds 和表格重新取数，不手抄结果。

## 数据与隐私

本仓库不包含：

- 原始 count 矩阵或 FPKM 表；
- 样本编号、分组表或测序数据；
- 参考基因组、GFF、eggNOG / KEGG 数据库；
- 差异基因表、富集结果、图片或研究结论。

运行时请通过环境变量或本地 `data/` 目录提供输入；这些路径已加入 `.gitignore`。

## 参考

- Love MI, Huber W, Anders S. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biology*. 2014.
- Zhu A, Ibrahim JG, Love MI. Heavy-tailed prior distributions for sequence count data: removing the noise and preserving large differences. *Bioinformatics*. 2019.
- Benjamini Y, Hochberg Y. Controlling the false discovery rate. *JRSS B*. 1995.
- Ritchie ME, et al. limma powers differential expression analyses for RNA-sequencing and microarray studies. *Nucleic Acids Research*. 2015.

