# ============================================================
# 09_summary.R  —  汇总报告 (关键统计表 + 候选靶标 + 可读报告)
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

## 2) 前置校验: 各中间 rds 必须存在, 否则明确指向上游步骤
load_rds <- function(name) {
  p <- file.path(rdsdir, name)
  if (!file.exists(p))
    stop(sprintf("缺失中间结果: %s\n  请先运行上游步骤生成该文件", p))
  readRDS(p)
}
f  <- load_rds("02_filtered.rds")
d  <- load_rds("04_deseq.rds")
tg <- load_rds("05_targets.rds")
en <- load_rds("07_enrichment.rds")

## 关键字段校验 (任一缺失都给出明确报错)
req_field <- function(obj, field, where)
  if (is.null(obj[[field]]))
    stop(sprintf("%s 缺少必填字段 '%s'", where, field))
req_field(f,  "count_mat", "02_filtered.rds")
req_field(d,  "n_OE",      "04_deseq.rds")
req_field(d,  "n_RNAi",    "04_deseq.rds")
req_field(d,  "n_OR",      "04_deseq.rds")
req_field(tg, "act",       "05_targets.rds")
req_field(tg, "rep",       "05_targets.rds")
req_field(tg, "cand_TF","05_targets.rds")
req_field(en, "enr",       "07_enrichment.rds")

## 从 03 保存的参数读取阈值 (单一事实来源, 与 DEG 判定一致)
thr_p <- if (!is.null(d$params$padj_thresh)) d$params$padj_thresh else 0.05
thr_l <- if (!is.null(d$params$lfc_thresh))  d$params$lfc_thresh  else 1
def_deg <- sprintf("padj<%.3f & |log2FC|>=%.1f (收缩LFC)", thr_p, thr_l)

## 1) 富集显著项计数: 具名安全取值
##    - 表缺失 (key 不存在): 告警并记 NA, 避免静默误判
##    - 表存在但为空 (如靶标集为空导致无特征): 返回 0 (合法结果, 非缺失)
enr_count <- function(key, thr = 0.05) {
  e <- en$enr[[key]]
  if (is.null(e)) {
    warning(sprintf("富集表 %s 缺失, 显著数记为 NA (请检查 07_enrichment 输出)", key))
    return(NA_integer_)
  }
  if (nrow(e) == 0) return(0L)
  if (!"padj" %in% names(e)) stop(sprintf("富集表 %s 缺少 padj 列", key))
  sum(e$padj < thr, na.rm = TRUE)
}
len <- function(x) if (is.null(x)) NA_integer_ else length(x)

summary <- data.frame(
  Item = c("保留基因数", "OE vs CK DEG", "RNAi vs CK DEG", "OE vs RNAi DEG",
           "激活靶标 (OE↑ & RNAi↓)", "抑制靶标 (OE↓ & RNAi↑)",
           "激活靶标 (探索性宽松)", "抑制靶标 (探索性宽松)",
           "PFAM富集(激活)显著数", "PFAM富集(抑制)显著数",
           "TF家族富集(激活)显著数", "TF家族富集(抑制)显著数"),
  Value = c(nrow(f$count_mat), d$n_OE, d$n_RNAi, d$n_OR,
            len(tg$act), len(tg$rep), len(tg$relaxed_act), len(tg$relaxed_rep),
            enr_count("PFAM_activated"), enr_count("PFAM_repressed"),
            enr_count("TFfamily_activated"), enr_count("TFfamily_repressed")),
  `定义_阈值` = c(
    "过滤后保留基因 (count>=10 于 >=2 样本)", def_deg, def_deg, def_deg,
    "OE↑ 且 RNAi↓ 均显著 (双侧反向保守筛选)",
    "OE↓ 且 RNAi↑ 均显著 (双侧反向保守筛选)",
    "一侧显著 + 另一侧反向趋势 (探索性, 非主结论)",
    "一侧显著 + 另一侧反向趋势 (探索性, 非主结论)",
    "超几何检验 BH padj<0.05", "超几何检验 BH padj<0.05",
    "超几何检验 BH padj<0.05", "超几何检验 BH padj<0.05"),
  stringsAsFactors = FALSE)
write.csv(summary, file.path(outdir, "summary.csv"), row.names = FALSE)

## 4) 候选位点 Top3 独立导出 (便于补充材料引用)
##    05 的 cand_TF = OE_only 子集 (OE 显著变化 & RNAi 无显著变化), 不再是混合 score top-N
cand_top3 <- head(tg$cand_TF, 3)
for (cl in c("evidence", "log2FC_RNAi", "padj_RNAi", "baseMean"))   # 兼容旧版 05 输出
  if (is.null(cand_top3[[cl]])) cand_top3[[cl]] <- NA
write.csv(cand_top3, file.path(outdir, "tables", "candidate_TF_top3.csv"),
          row.names = FALSE)

## 5) 可读 Markdown 报告 (可直接放入附件/补充材料)
shrink_type <- if (!is.null(d$params$shrink_type)) d$params$shrink_type else "apeglm"

## GO 注释来源按 06/07 的实际产物动态描述
##   (原为写死的 "GO 来自 UniProt 同源转移注释", 与默认的目标物种原生注释不符)
go_rds <- file.path(rdsdir, "06_go.rds")
go_ann <- if (file.exists(go_rds)) readRDS(go_rds) else NULL
gs <- en$go_stats
go_src_txt <- if (!is.null(gs)) {
  src <- if (identical(gs$source, "local"))
           sprintf("本地注释文件 %s", if (is.null(gs$source_file)) "(路径未记录)" else gs$source_file)
         else "UniProt 同源转移注释 (联网抓取)"
  sprintf("  GO 来自 %s (%d 个基因 / %d 个 GO 术语); PFAM/TF 家族来自输入注释表。",
          src, gs$genes, gs$terms)
} else if (!is.null(go_ann) &&
           identical(if (is.null(go_ann$status)) "ok" else go_ann$status, "ok")) {
  sprintf("  GO 来自 %s; PFAM/TF 家族来自输入注释表。", go_ann$source)
} else {
  "  ⚠ 本次**未产出 GO 富集** (06 步骤未给出可用注释; 不等于“无富集”); PFAM/TF 家族来自输入注释表。"
}

md <- c(
  sprintf("# TF RNA-seq 分析汇总报告"),
  sprintf("_生成时间: %s_", Sys.time()),
  "",
  "## 1. 关键统计",
  "",
  "| 指标 | 值 | 定义 / 阈值 |",
  "|------|----|--------------|",
  sprintf("| %s | %s | %s |", summary$Item, summary$Value, summary[["定义_阈值"]]),
  "",
  "## 2. 候选位点 (Top 3, 按 score 排序; 仅 OE_only: OE 显著变化 & RNAi 无显著变化)",
  "",
  "完整候选 (含全部 score top-N、证据分级、双侧 log2FC/padj) 见 `tables/candidate_reciprocal_score.csv`。",
  "",
  if (nrow(cand_top3) > 0)
    c("| gene_id | gene_name | Family | log2FC(OE/CK) | padj(OE) | log2FC(RNAi/CK) | padj(RNAi) | score |",
      "|---------|-----------|--------|---------------|----------|-----------------|------------|-------|",
      sprintf("| %s | %s | %s | %.3f | %.3g | %.3f | %.3g | %.3f |",
              cand_top3$gene_id, cand_top3$gene_name, cand_top3$Family,
              cand_top3$log2FoldChange, cand_top3$padj,
              cand_top3$log2FC_RNAi, cand_top3$padj_RNAi, cand_top3$score))
  else "_(无候选)_",
  "",
  "## 3. 方法摘要",
  "",
  sprintf("- **差异表达**: DESeq2 (Love et al., 2014), 设计 `%s`;",
          if (!is.null(d$params$design)) d$params$design else "~ condition"),
  sprintf("  log2FC 经 lfcShrink (type=%s) 收缩以降低低表达基因偏差。", shrink_type),
  sprintf("- **DEG 阈值**: padj < %.3f 且 |log2FC| >= %.1f。", thr_p, thr_l),
  "- **靶标筛选**: 激活=OE显著上调且RNAi显著下调; 抑制=反向 (双侧显著反向 = 保守策略)。",
  "- **富集**: 超几何检验, 背景 = \"有资格被判为 DEG 的基因\"(逐分析取 padj 非 NA 者; 不能判 DEG 的基因不进背景), BH 校正 (Benjamini & Hochberg, 1995);",
  "  背景口径与敏感性对照见 `tables/enrichment_methods.txt` 与 `tables/enrich_background_sensitivity.csv`;",
  go_src_txt,
  "- **详见**: `tables/references.md`(总文献), `tables/preprocessing_methods.txt`, `tables/DEG_methods.txt`, `tables/targets_methods.txt`, `tables/viz_methods.txt`。",
  "",
  "_注: 激活/抑制靶标集为双侧显著反向的保守结果; 若为空, 见探索性宽松集与候选位点表。_")
writeLines(md, file.path(outdir, "summary_report.md"))

## 主文献清单 (tables/*_methods.txt 引用均回溯至此)
refs <- c(
  "# References - TF RNA-seq analysis pipeline",
  "",
  "Master list cited across tables/*_methods.txt and the summary report.",
  "",
  "## Software / tools",
  "- DESeq2 - Love MI, Huber W, Anders S. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biol. 2014;15(12):550. DOI:10.1186/s13059-014-0550-8 | PMID:25516281",
  "- apeglm (lfcShrink) - Zhu A, Ibrahim JG, Love MI. Heavy-tailed prior distributions for sequence count data. Genome Biol. 2019;20(1):39. DOI:10.1186/s13059-019-1659-x | PMID:30799866",
  "- limma (venn) - Ritchie ME, Phipson B, Wu D, et al. Nucleic Acids Res. 2015;43(7):e47. DOI:10.1093/nar/gkv007 | PMID:25605792",
  "- edgeR (filtering) - Chen Y, Lun ATL, Smyth GK. Stat Anal Genomics Data Protoc. 2016:51-74.",
  "- voom/limma (filtering) - Law CW, Chen Y, Shi W, Smyth GK. Genome Biol. 2014;15(2):R29. DOI:10.1186/gb-2014-15-2-r29 | PMID:24485249",
  "- ggplot2 - Wickham H. ggplot2: Elegant Graphics for Data Analysis. Springer; 2016.",
  "- pheatmap - Kolde R. pheatmap: Pretty Heatmaps. R package; 2019.",
  "- STAR (aligner, upstream) - Dobin A, et al. Bioinformatics. 2013;29(1):15-21. DOI:10.1093/bioinformatics/bts635 | PMID:23104886",
  "- HISAT2 (aligner, upstream) - Kim D, et al. Nat Methods. 2015;12(4):357-360. DOI:10.1038/nmeth.3317 | PMID:25751142",
  "- featureCounts (counter, upstream) - Liao Y, et al. Bioinformatics. 2014;30(7):923-930. DOI:10.1093/bioinformatics/btt656 | PMID:24227677",
  "- HTSeq (counter, upstream) - Anders S, Pyl PT, Huber W. Bioinformatics. 2015;31(2):166-169. DOI:10.1093/bioinformatics/btu638 | PMID:25260700",
  "",
  "## Methods / statistics",
  "- BH FDR control - Benjamini Y, Hochberg Y. J R Stat Soc B. 1995;57(1):289-300.",
  "- Independent filtering - Bourgon R, Gentleman R, Huber W. Genome Biol. 2010;11(8):R94. DOI:10.1186/gb-2010-11-8-r94 | PMID:20727891",
  "- Enrichment background - Wijesooriya K, et al. Urgent need for consistent standards in functional",
  "  enrichment analysis. Brief Bioinform. 2022;23(3):bbac066. DOI:10.1093/bib/bbac066 | PMID:35272374",
  "- Enrichment background bias - Timmons JA, Szkop KJ, Gallagher IJ. Multiple sources of bias confound",
  "  functional enrichment analysis of global -omics data. Genome Biol. 2015;16:186. DOI:10.1186/s13059-015-0761-7 | PMID:26346307",
  "- Clustering / heatmap scaling - Eisen MB, et al. Proc Natl Acad Sci U S A. 1998;95(25):14863-14868. DOI:10.1073/pnas.95.25.14863 | PMID:9843981",
  "",
  "## Biology (TF family nomenclature)",
  "- Plant transcription factor families - Jin J, et al. PlantTFDB 4.0: toward a central hub for",
  "  transcription factors and regulatory interactions in plants. Nucleic Acids Res.",
  "  2017;45(D1):D1040-D1045. DOI:10.1093/nar/gkw982 | PMID:27924042",
  "- NOTE: the per-gene TF family assignment (column 'Family' of the input annotation table) is",
  "  third-party/vendor annotation. Its generating tool and version are NOT recorded anywhere in",
  "  this project, so they must be stated in the Methods. PlantTFDB is cited for the family",
  "  nomenclature only (bHLH / NAC / MYB_related / M-type_MADS / Nin-like ...), not as the source",
  "  of these particular assignments.",
  "",
  "## Notes on custom (non-standard) steps",
  "- Target score = log2FC(OE/CK) - log2FC(RNAi/CK) is a project-specific heuristic, not a published statistic.",
  "- Reciprocal OE/RNAi screen is a standard genetic strategy for TF target identification; DE stats are standard (DESeq2 + BH).",
  "- GO annotation transferred by homology from UniProt accessions in gene descriptions; same hypergeometric + BH framework as PFAM/TF.")
writeLines(refs, file.path(outdir, "tables", "references.md"))

cat(sprintf("\n[%s] 分析完成!\n", Sys.time()))
cat("---- 汇总 ----\n"); print(summary)
cat(sprintf("\n候选 TF 位点(Top3, 基于 OE↑/RNAi↓ 差异最大):\n"))
print(cand_top3)
cat(sprintf("\n结果已保存至: %s\n  - summary.csv / summary_report.md\n  - tables/candidate_TF_top3.csv\n",
            outdir))
