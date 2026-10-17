# ============================================================
# 00_config.R  —  公共配置 (被所有步骤 source)
# 环境/库 + 路径 + 常量; 幂等, 可重复 source
# 设计要点: 自动安装缺失包 / 校验样本列序 / 路径可经环境变量覆盖
# ============================================================
## 用户可写库(放置 pheatmap/svglite 等); 置于搜索路径首位
user_lib <- Sys.getenv("R_LIBS_USER", "~/R/library")
dir.create(user_lib, showWarnings = FALSE, recursive = TRUE)
.libPaths(c(user_lib, .libPaths()))

## ---------- 1) 缺失包自动安装 (新环境可直接 source) ----------
ensure_pkgs <- function() {
  ## CRAN 包
  cran <- c("ggplot2", "pheatmap", "RColorBrewer", "svglite")
  for (p in cran)
    if (!requireNamespace(p, quietly = TRUE))
      install.packages(p, repos = "https://cloud.r-project.org",
                       lib = .libPaths()[1])
  ## Bioconductor 包 (DESeq2, limma, apeglm 用于 lfcShrink)
  if (!requireNamespace("DESeq2",   quietly = TRUE) ||
      !requireNamespace("limma",    quietly = TRUE) ||
      !requireNamespace("apeglm",   quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager", repos = "https://cloud.r-project.org",
                       lib = .libPaths()[1])
    BiocManager::install(c("DESeq2", "limma", "apeglm"),
                         ask = FALSE, update = FALSE)
  }
}
ensure_pkgs()

## 仅屏蔽包启动消息(保留后续运行期 warning, 便于发现版本/函数冲突)
suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
  library(limma)
})

## 固定随机种子, 保证任何随机步骤(抽样/置换)可复现
set.seed(123)

## 确保运行期 warning 立即显示 (不关闭报错/告警提示; stop() 错误本就中断流程)
##   注意: 上面 suppressPackageStartupMessages 仅屏蔽包加载消息, 不影响 warning/error
options(warn = 1)

## ---------- 2) 路径 (可用环境变量覆盖, 便于迁移) ----------
base   <- Sys.getenv("RNASEQ_BASE",   "~/rna-seq")
countf <- Sys.getenv("RNASEQ_COUNT", file.path(base, "gene_count.xls"))
outdir <- file.path(base, "results")
rdsdir <- file.path(outdir, "rds")

dir.create(outdir,        showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(outdir, "tables"),  showWarnings = FALSE)
dir.create(file.path(outdir, "figures"), showWarnings = FALSE)
dir.create(rdsdir,        showWarnings = FALSE)

## ---------- 3) 样本与分组 (顺序须与 count 矩阵列序一致) ----------
## 原始数据列序: KD1,KD2,KD3, OE1,OE2,OE3, CK1,CK2
## 注意: CK 仅 2 个生物学重复 -> DESeq2 离散度估计偏不稳, 差异基因假阳性偏高,
##       结果(尤其 CK 相关对比)解读需谨慎; 此为实验设计层面局限, 非代码可修。
sample_cols <- c("KD1","KD2","KD3","OE1","OE2","OE3","CK1","CK2")
cond <- factor(c(rep("RNAi", 3), rep("OE", 3), rep("CK", 2)),
               levels = c("CK", "RNAi", "OE"))

## 样本列校验: 只需确保 gene_id 与全部 sample_cols 列名存在。
## 注意顺序无关: 下游 df[, sample_cols] 会按 sample_cols 顺序重排,
## 因此文件列序变化不影响分组, 无需强制顺序一致 (避免误杀合法文件)。
validate_samples <- function() {
  if (!file.exists(countf)) stop("count 文件不存在: ", countf)
  hdr <- strsplit(readLines(countf, n = 1, warn = FALSE), "\t")[[1]]
  if (!("gene_id" %in% hdr)) stop("count 文件缺少 'gene_id' 列")
  miss <- setdiff(sample_cols, hdr)
  if (length(miss) > 0)
    stop("count 文件缺少样本列: ", paste(miss, collapse = ", "),
         "\n  现有列: ", paste(hdr, collapse = ", "))
  invisible(TRUE)
}
validate_samples()

cat(sprintf("[%s] 配置加载完成 (outdir=%s)\n", Sys.time(), outdir))
