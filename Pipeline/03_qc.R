# ============================================================
# 03_qc.R  —  样本质控: PCA (置于 04 差异分析之前)
#   先看样本结构再解读 DE。对 DESeq2 vst(blind=TRUE) 后的表达矩阵,
#   取方差最高 N 个基因做 PCA。blind=TRUE 不使用设计公式, 展示样本间
#   总变异(含潜在批次效应), 是 RNA-seq 样本 QC 的标准做法。
#   文献: Love, Huber & Anders 2014 (DESeq2 / vst)。方法详见 tables/viz_methods.txt。
# ============================================================
script_dir <- (function() {
  args <- commandArgs(trailingOnly = FALSE)
  fa <- grep("^--file=", args, value = TRUE)
  if (length(fa)) return(dirname(sub("^--file=", "", fa[1])))
  for (i in sys.nframe():1) {
    of <- tryCatch(sys.frame(i)$ofile, error = function(e) NULL)
    if (!is.null(of) && of != "") return(dirname(of))
  }
  "."
})()
source(file.path(script_dir, "00_config.R"))

f <- readRDS(file.path(rdsdir, "02_filtered.rds"))
count_mat <- f$count_mat
pca_ntop <- as.integer(Sys.getenv("PCA_NTOP", "500"))

## 仅为 QC 构建 dds (blind=TRUE 忽略设计; 与 04_deseq 的 dds 相互独立)
coldata <- data.frame(condition = cond, row.names = sample_cols)
dds <- DESeqDataSetFromMatrix(countData = count_mat, colData = coldata, design = ~ condition)
vsd <- vst(dds, blind = TRUE)
mat <- assay(vsd)
rv  <- rowSums((mat - rowMeans(mat))^2) / (ncol(mat) - 1)
top <- head(order(rv, decreasing = TRUE), pca_ntop)
pca <- prcomp(t(mat[top, ]), scale. = FALSE)
pcaData <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2],
                      condition = cond, name = sample_cols)
pct <- round(pca$sdev^2 / sum(pca$sdev^2) * 100, 1)
p_pca <- ggplot(pcaData, aes(x = PC1, y = PC2, color = condition, label = name)) +
  geom_point(size = 3) + geom_text(vjust = -1.2, size = 3) +
  xlab(sprintf("PC1: %s%% variance", pct[1])) +
  ylab(sprintf("PC2: %s%% variance", pct[2])) +
  ggtitle("PCA (vst, sample QC before DE)") +
  theme_bw(base_family = "Arial", base_size = 12) +
  theme(axis.text = element_text(size = 11),
        axis.title = element_text(size = 12),
        plot.title = element_text(face = "bold", hjust = 0.5, size = 13)) +
  scale_color_manual(values = c(CK = "grey", RNAi = "blue", OE = "red"))
ggsave(file.path(outdir, "figures", "PCA.svg"), p_pca, width = 6.5, height = 5.5)
## PCA.svg 字体统一为 Arial (cairo 本地用 Liberation Sans 度量, 与 Arial 一致)
{ fp <- file.path(outdir, "figures", "PCA.svg")
  t <- readLines(fp, warn = FALSE)
  writeLines(gsub('font-family: "[^"]*"', 'font-family: "Arial"', t), fp) }
write.csv(pcaData, file.path(outdir, "tables", "PCA_coords.csv"), row.names = FALSE)
saveRDS(list(pcaData = pcaData, pct = pct), file.path(rdsdir, "03_qc.rds"))

cat(sprintf("[%s] QC 完成: PCA PC1=%.1f%% / PC2=%.1f%% (top %d 变异基因)。\n  请先查看 figures/PCA.svg 的样本聚类/批次结构, 再进行 04 差异分析。\n",
            Sys.time(), pct[1], pct[2], pca_ntop))
