# ============================================================
# 08_visualization.R  —  火山 / 维恩 / 热图 / 富集条形图 / 靶标散点
#   (PCA / 样本质控已在 03_qc.R 完成, 本步只画 DE 相关图)
# 注意标准化用途:
#   * 热图用 DESeq2 归一化 count 的 log2(count+1) -> 仅展示表达量, 非 vst
#   * PCA(在 03_qc)才用 vst 方差稳定转换
# ============================================================
## 本脚本所在目录 (确保无论从何处调用都能找到同目录文件)
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

d  <- readRDS(file.path(rdsdir, "04_deseq.rds"))
en <- readRDS(file.path(rdsdir, "07_enrichment.rds"))
dds <- d$dds
## 显著性阈值取自 04 的 params (单一事实来源), 与 DEG 判定保持一致
thr_p <- if (!is.null(d$params$padj_thresh)) d$params$padj_thresh else 0.05
thr_l <- if (!is.null(d$params$lfc_thresh))  d$params$lfc_thresh  else 1

## 可视化参数 (可经环境变量覆盖, 调参不改源码)
heat_top  <- as.integer(Sys.getenv("HEATMAP_TOP", "100"))
spec_top  <- as.integer(Sys.getenv("SPEC_TOP",     "12"))
go_top    <- as.integer(Sys.getenv("GO_TOP",       "10"))
kegg_top  <- as.integer(Sys.getenv("KEGG_TOP",     "15"))
plots     <- list()   # 缓存 ggplot 对象, 末尾落盘供后续微调

norm_counts <- log2(counts(dds, normalized = TRUE) + 1)

## (PCA / 样物质控已移至 03_qc.R, 在差异分析前完成; 本步只画 DE 相关图)

## 7.1 火山图 (三组对比)
volcano <- function(tab, title, fname) {
  ## 先剔除 padj/log2FC 缺失行 (独立过滤剔除的基因), 避免绘图时产生 NA 告警
  tab <- tab[!is.na(tab$padj) & !is.na(tab$log2FoldChange), ]
  tab$signif <- (!is.na(tab$padj) & tab$padj < thr_p &
                 !is.na(tab$log2FoldChange) & abs(tab$log2FoldChange) >= thr_l)
  ## 坐标轴自适应: 以数据最大 |log2FC| 为界(下限 10), 超出则告警避免丢失极端差异点
  xmax <- max(10, ceiling(max(abs(tab$log2FoldChange), na.rm = TRUE)))
  nclip <- sum(abs(tab$log2FoldChange) > xmax, na.rm = TRUE)
  if (nclip > 0)
    warning(sprintf("%s: %d 个基因 |log2FC| > %g 已贴边显示 (未被丢弃)", title, nclip, xmax))
  p <- ggplot(tab, aes(x = log2FoldChange, y = -log10(padj), color = signif)) +
    geom_point(size = 0.6, alpha = 0.5) +
    scale_color_manual(values = c("grey", "red")) +
    geom_hline(yintercept = -log10(thr_p), linetype = "dashed") +
    geom_vline(xintercept = c(-thr_l, thr_l), linetype = "dashed") +
    ggtitle(title) + theme_bw(base_family = "Arial") + xlim(-xmax, xmax)
  ggsave(file.path(outdir, "figures", fname), p, width = 6, height = 5)
  p
}
plots[["Volcano_OE_vs_CK"]]   <- volcano(d$t_OE_CK,   "OE vs CK",   "Volcano_OE_vs_CK.svg")
plots[["Volcano_RNAi_vs_CK"]] <- volcano(d$t_RNAi_CK, "RNAi vs CK", "Volcano_RNAi_vs_CK.svg")
plots[["Volcano_OE_vs_RNAi"]] <- volcano(d$t_OE_RNAi, "OE vs RNAi", "Volcano_OE_vs_RNAi.svg")

## 7.3 维恩图 (DEG 重叠) — 用 vennCounts 构造标准结构 (原脚本直接传 3x3 矩阵, 语义错误)
## 注意: limma::vennCounts 要求 行=元素(基因), 列=集合(对比)
sets <- list("OE_vs_CK" = d$deg_OE, "RNAi_vs_CK" = d$deg_RNAi, "OE_vs_RNAi" = d$deg_OR)
## 防御: 各 DEG 集合及并集不应含重复 gene_id (重复会令 0/1 矩阵计数错乱)
for (s in names(sets))
  if (any(duplicated(sets[[s]])))
    stop(sprintf("维恩图: %s 含重复 gene_id, 请检查 04 输出", s))
all_ids <- union(union(d$deg_OE, d$deg_RNAi), d$deg_OR)
mem <- matrix(0L, nrow = length(all_ids), ncol = 3,
              dimnames = list(NULL, names(sets)))
for (i in seq_along(sets)) mem[which(all_ids %in% sets[[i]]), i] <- 1
vc <- limma::vennCounts(mem)
svg(file.path(outdir, "figures", "Venn_DEG.svg"), width = 6, height = 6)
limma::vennDiagram(vc, names = names(sets), circle.col = c("red", "blue", "green"))
dev.off()

## 7.4 热图 (Top DEGs)
top_ids <- en$union_deg
if (is.null(top_ids) || length(top_ids) == 0) {
  warning("热图: en$union_deg 为空, 跳过 Top DEG 热图")
} else {
  if (length(top_ids) > heat_top) {
    ord <- order(-abs(d$t_OE_CK$log2FoldChange[match(en$union_deg, d$t_OE_CK$gene_id)]))
    top_ids <- en$union_deg[ord][seq_len(heat_top)]
  }
  mat <- norm_counts[top_ids, , drop = FALSE]
  ## 逐基因 z-score 标准化 (减行均值, 除以行标准差; 标准差为0/NA 则置1避免 NaN)
  ## —— 表达量热图的常规逐基因归一化 (同 pheatmap scale="row"), 便于跨基因比较表达模式
  rmean <- rowMeans(mat); rsd <- apply(mat, 1, sd)
  rsd[rsd == 0 | is.na(rsd)] <- 1
  mat <- (mat - rmean) / rsd
  if (length(top_ids) > 1) {
    svg(file.path(outdir, "figures", "Heatmap_topDEG.svg"), width = 11, height = 10)
    ph <- pheatmap(mat, scale = "none", show_rownames = TRUE, fontsize_row = 5,
                   annotation_col = data.frame(Condition = cond, row.names = colnames(mat)),
                   main = sprintf("Top %d DEGs (row z-score)", length(top_ids)),
                   fontsize_col = 10)
    dev.off()
    ## 导出"图中行序(上->下) -> 基因"映射 + OE/RNAi 统计与注释, 便于定位每一行是哪个基因
    ord_ids <- if (!is.null(ph$tree_row)) ph$tree_row$labels[ph$tree_row$order] else rownames(mat)
    hm <- data.frame(row_position = seq_along(ord_ids), gene_id = ord_ids, stringsAsFactors = FALSE)
    toe <- d$t_OE_CK; tri <- d$t_RNAi_CK
    hm$log2FC_OE   <- round(toe$log2FoldChange[match(hm$gene_id, toe$gene_id)], 3)
    hm$padj_OE     <- signif(toe$padj[match(hm$gene_id, toe$gene_id)], 3)
    hm$log2FC_RNAi <- round(tri$log2FoldChange[match(hm$gene_id, tri$gene_id)], 3)
    hm$padj_RNAi   <- signif(tri$padj[match(hm$gene_id, tri$gene_id)], 3)
    anno <- readRDS(file.path(rdsdir, "02_filtered.rds"))$anno   # 注释(Family/描述)
    hm <- merge(hm, anno[, intersect(c("gene_id","Family","gene_description"), colnames(anno))],
                by = "gene_id", all.x = TRUE)
    hm <- hm[order(hm$row_position), ]
    write.csv(hm, file.path(outdir, "tables", "heatmap_gene_order.csv"), row.names = FALSE)
  }
}

## 7.5 富集条形图 — 激活 vs 抑制 特异性 (科学层面: 2x2 Fisher, 稀有域已过滤)
##   橙=激活富集, 蓝=抑制富集; 标注 in_activated/in_repressed; 红虚线=p<0.05
##   取 pvalue 最小 Top12 候选; 若均不显著, 条全部低于红线 -> 如实呈现"无特异性富集"
for (nm in c("PFAM_spec", "TFfamily_spec")) {
  e <- en$enr[[nm]]
  if (is.null(e) || nrow(e) == 0) next
  ## 防御: 特异性富集表必须含标注所需列, 否则报错而非静默 NA
  if (!all(c("in_activated", "in_repressed", "direction") %in% names(e)))
    stop(sprintf("富集条形图 %s 缺少 in_activated/in_repressed/direction 列", nm))
  e <- head(e[order(e$pvalue), ], spec_top)
  e$feature   <- factor(e$feature, levels = rev(e$feature))
  e$bar_color <- ifelse(e$direction == "activated", "orange",
                  ifelse(e$direction == "repressed", "blue", "grey50"))
  ymax <- max(-log10(e$pvalue)) * 1.18
  npass <- sum(e$padj < 0.05, na.rm = TRUE)
  p <- ggplot(e, aes(x = feature, y = -log10(pvalue), fill = bar_color)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = sprintf("%d/%d", in_activated, in_repressed)),
              hjust = -0.1, size = 3, color = "black") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
    scale_fill_identity() +
    coord_flip() + expand_limits(y = c(0, ymax)) +
    xlab("domain (in_activated/in_repressed)") +
    ggtitle(sprintf("Specificity candidates (act vs rep): %s  [BH padj<0.05: %d]",
                    nm, npass)) + theme_bw()
  ggsave(file.path(outdir, "figures", sprintf("Bar_%s.svg", nm)), p,
         width = 7, height = 5)
  plots[[sprintf("Bar_%s", nm)]] <- p
}

## 7.5 GO 富集气泡图 (出版级; 参照 clusterProfiler/ggplot2 bubble 风格)
##   x = GeneRatio(富集比例), y = GO term, 点大小 = Count(基因数), 颜色 = -log10(padj) 绿->红
##   目标物种原生注释较密集, 按 raw p 取 TopN 展示"最富集的 term"; 是否 BH 显著在标题注明
wrap_term <- function(s, w = 40)
  sapply(as.character(s), function(x) paste(strwrap(x, width = w), collapse = "\n"))
go_bubble <- function(e, title, fname, topn = 12) {
  if (is.null(e) || nrow(e) == 0) return(invisible(NULL))
  ## 去掉已废弃 GO 术语; 用 is.na() 兜底 —— 否则 NA 作为索引会插入整行 NA
  e <- e[is.na(e$go_term) | !grepl("^obsolete", e$go_term, ignore.case = TRUE), ]
  ## 按本体(BP/CC/MF)分组: 每类各取 topn, 分面展示, 使三类各自成 panel
  e$ns <- factor(e$namespace, levels = c("BP","CC","MF"),
                 labels = c("Biological process","Cellular component","Molecular function"))
  e <- e[order(e$pvalue), ]
  e <- do.call(rbind, lapply(split(e, e$ns), function(x) head(x, topn)))
  e <- e[order(e$ns, e$pvalue), ]
  e$GeneRatio <- e$in_test / e$test_total
  e$neglog10 <- -log10(pmax(e$padj, .Machine$double.xmin))
  wrapped <- wrap_term(e$go_term)
  e$lab <- factor(wrapped, levels = rev(unique(wrapped)))   # 每个分面内最显著在上
  npass <- sum(e$padj < 0.05, na.rm = TRUE)
  p <- ggplot(e, aes(x = GeneRatio, y = lab)) +
    geom_point(aes(size = in_test, colour = neglog10), alpha = 0.9) +
    scale_colour_gradient(low = "#2CB42C", high = "#D62728", name = "-log10(padj)") +
    scale_size_continuous(range = c(2.5, 6.5), name = "Count") +
    labs(x = "Gene ratio", y = NULL,
         title = sprintf("%s   (top %d per ontology; BH padj<0.05: %d)", title, topn, npass)) +
    facet_grid(rows = vars(ns), scales = "free_y", space = "free_y") +
    theme_bw(base_family = "Arial", base_size = 12) +
    theme(axis.text.y = element_text(size = 9),
          axis.text.x = element_text(size = 11),
          plot.title = element_text(face = "bold", hjust = 0.5, size = 13),
          strip.text = element_text(face = "bold", size = 11),
          strip.background = element_rect(fill = "grey92", colour = NA),
          panel.spacing = unit(0.4, "lines"),
          panel.grid.major = element_line(linewidth = 0.2, colour = "grey88"),
          panel.grid.minor = element_line(linewidth = 0.1, colour = "grey92"),
          legend.position = "right")
  ggsave(file.path(outdir, "figures", fname), p,
         width = 9, height = max(4.5, nrow(e) * 0.30 + 1.2))
  p
}
for (nm in names(en$enr)[grepl("^GO_", names(en$enr))])
  plots[[sprintf("Bubble_%s", nm)]] <-
    go_bubble(en$enr[[nm]], gsub("^GO_", "GO: ", nm), sprintf("Bubble_%s.svg", nm), go_top)

## KEGG 通路气泡图(同 GO 风格; 单面板, 无本体分层)
for (nm in names(en$enr)[grepl("^KEGG_", names(en$enr))]) {
  e <- en$enr[[nm]]
  if (is.null(e) || nrow(e) == 0) next
  e$nm_txt <- ifelse(is.na(e$name) | !nzchar(e$name), e$feature, e$name)
  e <- e[order(e$pvalue), ]
  e <- head(e, kegg_top)
  e$GeneRatio <- e$in_test / e$test_total
  e$neglog10 <- -log10(pmax(e$padj, .Machine$double.xmin))
  e$lab <- factor(wrap_term(paste0(e$nm_txt, " [", e$feature, "]")),
                  levels = rev(wrap_term(paste0(e$nm_txt, " [", e$feature, "]"))))
  npass <- sum(e$padj < 0.05, na.rm = TRUE)
  p <- ggplot(e, aes(x = GeneRatio, y = lab)) +
    geom_point(aes(size = in_test, colour = neglog10), alpha = 0.9) +
    scale_colour_gradient(low = "#2CB42C", high = "#D62728", name = "-log10(padj)") +
    scale_size_continuous(range = c(2.5, 6.5), name = "Count") +
    labs(x = "Gene ratio", y = NULL,
         title = sprintf("KEGG: %s   (top %d; BH padj<0.05: %d)", nm, kegg_top, npass)) +
    theme_bw(base_family = "Arial", base_size = 12) +
    theme(axis.text.y = element_text(size = 9), axis.text.x = element_text(size = 11),
          plot.title = element_text(face = "bold", hjust = 0.5, size = 13),
          panel.grid.major = element_line(linewidth = 0.2, colour = "grey88"),
          panel.grid.minor = element_line(linewidth = 0.1, colour = "grey92"),
          legend.position = "right")
  ggsave(file.path(outdir, "figures", sprintf("Bubble_%s.svg", nm)), p,
         width = 9, height = max(4.5, nrow(e) * 0.30 + 1.2))
  plots[[sprintf("Bubble_%s", nm)]] <- p
}

## 7.6 靶标方向散点图 (OE vs CK 与 RNAi vs CK 的反向一致性)
##   激活靶标: OE上调(x>0) 且 RNAi下调(y<0) -> 右下象限 (橙)
##   抑制靶标: OE下调(x<0) 且 RNAi上调(y>0) -> 左上象限 (蓝)
tg <- readRDS(file.path(rdsdir, "05_targets.rds"))
oe <- d$t_OE_CK; ri <- d$t_RNAi_CK
## 防御: 两张差异表的 gene_id 必须唯一, 否则 match 会引入 NA/错位
if (any(duplicated(oe$gene_id)) || any(duplicated(ri$gene_id)))
  stop("靶标方向散点: t_OE_CK / t_RNAi_CK 含重复 gene_id, 请检查 04/05 输出")
m_ri <- match(oe$gene_id, ri$gene_id)
if (any(is.na(m_ri)))
  warning("靶标方向散点: 部分 OE 基因在 RNAi 表中缺失, 对应点已忽略")
sdf <- data.frame(
  gene_id = oe$gene_id,
  oe_lfc  = oe$log2FoldChange,
  ri_lfc  = ri$log2FoldChange[m_ri]
)
sdf$class <- "other"
sdf$class[sdf$gene_id %in% tg$act] <- "activated"
sdf$class[sdf$gene_id %in% tg$rep] <- "repressed"
sdf <- sdf[!is.na(sdf$oe_lfc) & !is.na(sdf$ri_lfc), ]
sdf$class <- factor(sdf$class, levels = c("other", "activated", "repressed"))

## 坐标上限自适应: 固定 ±10 会把 |LFC|>10 的点(含最高分候选)裁到画布外, 图标会整块消失
lim <- max(10, ceiling(max(abs(c(sdf$oe_lfc, sdf$ri_lfc)), na.rm = TRUE)))
nclip <- sum(abs(c(sdf$oe_lfc, sdf$ri_lfc)) > lim, na.rm = TRUE)
if (nclip > 0)
  warning(sprintf("靶标散点: %d 个基因 |log2FC| > %g 已贴边显示 (未被丢弃)", nclip, lim))

## 候选位点 (05 的 OE_only 子集, 见 tables/candidate_TF_locus.csv): 图注动态生成,
## 不再把 "其 OE padj=NA" 之类的单次结果写死在标题里
cand_top <- if (!is.null(tg$cand_TF) && nrow(tg$cand_TF) > 0) tg$cand_TF[1, ] else NULL
cg <- if (is.null(cand_top)) NA_character_ else cand_top$gene_id
cand_subtitle <- if (!is.null(cand_top) && cg %in% sdf$gene_id) {
  ev <- if (is.null(cand_top$evidence))    "n/a"     else cand_top$evidence
  rl <- if (is.null(cand_top$log2FC_RNAi)) NA_real_  else cand_top$log2FC_RNAi
  rp <- if (is.null(cand_top$padj_RNAi))   NA_real_  else cand_top$padj_RNAi
  sprintf("红菱形: score 最高候选 %s [%s] (OE log2FC=%.2f, padj=%.2g | RNAi log2FC=%.2f, padj=%.2g)",
          cg, ev, cand_top$log2FoldChange, cand_top$padj, rl, rp)
} else "无 OE_only 候选位点 (candidate_TF_locus.csv 为空)"

p_dir <- ggplot(sdf, aes(x = oe_lfc, y = ri_lfc, color = class)) +
  geom_point(size = 0.5, alpha = 0.6) +
  scale_color_manual(values = c(other = "grey70",
                                activated = "orange",
                                repressed = "blue")) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  coord_cartesian(xlim = c(-lim, lim), ylim = c(-lim, lim)) +
  xlab("log2FC (OE vs CK)") + ylab("log2FC (RNAi vs CK)") +
  ggtitle("Target direction: OE vs CK vs RNAi vs CK") +
  labs(subtitle = cand_subtitle) +
  theme_bw(base_family = "Arial")
## 标注候选位点 (05 已保证其为 OE_only: OE 显著变化而 RNAi 无显著变化)
if (!is.null(cand_top)) {
  cp  <- sdf[sdf$gene_id == cg, ]
  cp$lab <- cg        # 标签写进数据本身: 否则缓存的 ggplot 对象依赖脚本变量 cg,
                      # 之后在别的会话里重绘(如 Report/03_figures.R)会报 object 'cg' not found
  if (nrow(cp) > 0) {
    p_dir <- p_dir +
      geom_point(data = cp, aes(x = oe_lfc, y = ri_lfc),
                 color = "red", size = 3, shape = 18) +
      geom_text(data = cp, aes(x = oe_lfc, y = ri_lfc, label = lab),
                vjust = -1, color = "red", size = 3)
  }
}
ggsave(file.path(outdir, "figures", "Scatter_target_direction.svg"), p_dir,
       width = 6, height = 6)
plots[["Scatter_target_direction"]] <- p_dir

## 缓存 ggplot 对象 (PCA/火山/富集/散点), 便于后续微调图样而不必重跑全流程
saveRDS(plots, file.path(rdsdir, "08_plots.rds"))

## 落盘可视化方法说明 (含文献), 保证文档随流水线再生、不被丢失
## PCA 方差贡献从 03_qc.rds 现读 (原为写死字符串, 换样本/换阈值后不会同步)
qc_rds <- file.path(rdsdir, "03_qc.rds")
qc <- if (file.exists(qc_rds)) readRDS(qc_rds) else NULL
pc <- if (!is.null(qc) && length(qc$pct) >= 3) qc$pct[1:3] else rep(NA_real_, 3)
pca_line <- if (!any(is.na(pc))) {
  sprintf("   - PC1=%.1f%%, PC2=%.1f%%, PC3=%.1f%% for this dataset; coords in tables/PCA_coords.csv.",
          pc[1], pc[2], pc[3])
} else {
  "   - See tables/PCA_coords.csv for this dataset's variance contributions (03_qc.rds unavailable)."
}
go_ns_note <- if (!is.null(en$go_stats) && !is.null(en$go_stats$terms_no_ns) && en$go_stats$terms_no_ns > 0) {
  sprintf("   - NOTE: %d GO terms of this annotation carry no BP/CC/MF namespace and are NOT plotted (they are retained in tables/enrich_GO_*.csv).",
          en$go_stats$terms_no_ns)
} else {
  "   - NOTE: only GO terms carrying a BP/CC/MF namespace are plotted."
}

viz_methods <- c(
  "Visualization methods (TF RNA-seq)",
  "=====================================",
  "All plots generated with ggplot2 (Wickham 2016) and pheatmap (Kolde 2019);",
  "DESeq2 normalization and variance-stabilizing transformation via DESeq2 (Love, Huber & Anders 2014).",
  "Significance thresholds (padj / |log2FC|) are read from 04_deseq.rds$params to stay consistent",
  "with the differential-expression calling.",
  "",
  "1) PCA",
  "   - vst (variance-stabilizing transformation, blind=TRUE) before PCA [Love et al. 2014];",
  "     prcomp() on top-N highest-variance genes (default N=500, env PCA_NTOP), scale.=FALSE.",
  pca_line,
  "   Ref: Love MI, Huber W, Anders S. Genome Biol. 2014;15(12):550. DOI:10.1186/s13059-014-0550-8 | PMID:25516281",
  "",
  "2) Volcano plots",
  "   - x = log2 fold change, y = -log10(padj); threshold padj < thr_p AND |log2FC| >= thr_l (from 04 params).",
  "   - x-axis auto-expands to max |log2FC| (floor 10); warns if points exceed it (no silent clipping).",
  "",
  "3) Venn diagram (DEG overlap)",
  "   - limma::vennDiagram over the three significant DEG sets [Ritchie et al. 2015].",
  "   Ref: Ritchie ME, et al. Nucleic Acids Res. 2015;43(7):e47. DOI:10.1093/nar/gkv007 | PMID:25605792",
  "",
  "4) Top-DEG heatmap",
  "   - Genes: union of the three DEG sets, ranked by |log2FC(OE vs CK)| descending, top N",
  "     (N = HEATMAP_TOP, default 100). DESeq2 normalized counts, log2(count+1);",
  "     per-gene row z-score = (x-rowMean)/rowSD",
  "     (conventional per-gene scaling, == pheatmap scale='row'; Eisen et al. 1998).",
  "   NOTE: log2(normalized+1) is a DISPLAY transform only; NOT a vst; PCA uses vst.",
  "   Ref: Eisen MB, et al. Proc Natl Acad Sci U S A. 1998;95(25):14863-14868. DOI:10.1073/pnas.95.25.14863 | PMID:9843981",
  "",
  "5) Enrichment plots",
  "   - PFAM/TF/KEGG: hypergeometric, BH-adjusted; specificity (act vs rep) 2x2 Fisher.",
  "   - Universe = genes eligible for DEG calling (non-NA padj after DESeq2 filtering); genes that",
  "     cannot be called DE are excluded so that the hypergeometric sampling model holds. Background",
  "     choice and a sensitivity comparison are documented in tables/enrichment_methods.txt and",
  "     tables/enrich_background_sensitivity.csv.",
  "   - GO: bubble plot, x = GeneRatio (in_test/test_total), dot size = Count (in_test),",
  "     colour = -log10(padj) green->red; faceted by ontology (BP / CC / MF), top-N per ontology",
  "     by raw p; BH<0.05 count in title. (clusterProfiler/ggplot2 publication style; no target species",
  "     OrgDb -> custom hypergeometric.)",
  go_ns_note,
  "   - Rare features (bg<5) filtered; BH-adjusted [Benjamini & Hochberg 1995].",
  "   Ref: Benjamini Y, Hochberg Y. J R Stat Soc B. 1995;57(1):289-300.",
  "",
  "7) Fonts",
  "   - Every ggplot figure is built with base_family='Arial' (so the text metrics are taken from a",
  "     metric-compatible substitute, Liberation Sans), and the SVG font-family declarations are",
  "     rewritten to Arial at the end so journals/Illustrator render Arial. Heatmap (pheatmap) and",
  "     Venn (limma) are emitted by the cairo SVG device as vector glyph outlines (no font-family",
  "     declaration, i.e. font-independent).",
  "",
  "6) Target-direction scatter (OE vs CK vs RNAi vs CK)",
  "   - x=log2FC(OE/CK), y=log2FC(RNAi/CK); activated=lower-right, repressed=upper-left.",
  "   - Axis limits expand to the data range (floor +-10) so that extreme points are not clipped.",
  "   - The top-ranked OE_only candidate (05: significant change in OE, no significant change in",
  "     RNAi; tables/candidate_TF_locus.csv) is marked with a red diamond and labelled with its",
  "     gene id; its log2FC/padj on both sides are printed in the subtitle. Such a locus is an",
  "     over-expression/insertion candidate, NOT a reciprocal target; identity needs BLAST.",
  "",
  "Software: ggplot2 (Wickham 2016); pheatmap (Kolde 2019).")
writeLines(viz_methods, file.path(outdir, "tables", "viz_methods.txt"))

## 统一字体: 所有 SVG 的 font-family 声明改为 Arial (放在最末, 确保散点等后生成的图也覆盖)
##   本地 cairo 用 Liberation Sans 度量(与 Arial 度量一致); 改写声明便于期刊/Illustrator 按 Arial 渲染。
##   注: 热图(pheatmap)/维恩(limma)为路径化字形, 无 font-family, 本身不依赖字体、任何查看器一致。
force_arial_svg <- function() {
  svgs <- list.files(file.path(outdir, "figures"), pattern = "\\.svg$", full.names = TRUE)
  for (f in svgs) {
    txt <- readLines(f, warn = FALSE)
    txt <- gsub('font-family: "[^"]*"', 'font-family: "Arial"', txt)
    writeLines(txt, f)
  }
  length(svgs)
}
cat(sprintf("[%s] 已将 %d 个 SVG 字体统一为 Arial\n", Sys.time(), force_arial_svg()))

cat(sprintf("[%s] 可视化完成 (%d 张 ggplot 对象已缓存)\n", Sys.time(), length(plots)))
