# ============================================================
# 10_report.R  —  生成详细分析报告 ANALYSIS_REPORT.md (由流水线再生, 免受工作区重置删除)
#   所有数字均从产物现算 (tables/filter_stats.csv / rds/03_qc.rds / rds/04,05,06,07),
#   不再写死: 换阈值、换注释来源或加样本后, 报告数字会跟着变, 不会静默输出旧结论。
#   用法: Rscript Pipeline/10_report.R      (需先跑过 01 → 07)
#   注: 全流程不存在一键包装脚本, 各步骤单独运行 (见报告第 10 节)
# ============================================================
script_dir <- (function() {
  args <- commandArgs(trailingOnly = FALSE)
  fa <- grep("^--file=", args, value = TRUE)
  if (length(fa)) return(dirname(sub("^--file=", "", fa[1])))
  "."
})()
source(file.path(script_dir, "00_config.R"))

f  <- readRDS(file.path(rdsdir, "02_filtered.rds"))
d  <- readRDS(file.path(rdsdir, "04_deseq.rds"))
tg <- readRDS(file.path(rdsdir, "05_targets.rds"))
n_gene <- nrow(f$count_mat)
n_oe <- d$n_OE; n_ri <- d$n_RNAi; n_or <- d$n_OR
n_act <- length(tg$act); n_rep <- length(tg$rep)
n_ract <- length(tg$relaxed_act); n_rrep <- length(tg$relaxed_rep)
thr_p <- if (!is.null(d$params$padj_thresh)) d$params$padj_thresh else 0.05
thr_l <- if (!is.null(d$params$lfc_thresh))  d$params$lfc_thresh  else 1
shrink_type <- if (!is.null(d$params$shrink_type)) d$params$shrink_type else "apeglm"

## ============================================================
## 一、现算所有报告数字
## ============================================================
## 1) 过滤统计
fs <- read.csv(file.path(outdir, "tables", "filter_stats.csv"), stringsAsFactors = FALSE)
fk <- function(k) fs$value[match(k, fs$parameter)]
n_raw <- fk("n_before"); f_min_count <- fk("min_count"); f_min_samples <- fk("min_samples")

## 2) QC (03): 方差贡献 + 关键样本坐标
qc_rds <- file.path(rdsdir, "03_qc.rds")
qc <- if (file.exists(qc_rds)) readRDS(qc_rds) else NULL
pc <- if (!is.null(qc) && length(qc$pct) >= 3) qc$pct[1:3] else rep(NA_real_, 3)
pct_txt <- function(i) if (is.na(pc[i])) "未计算" else sprintf("%.1f%%", pc[i])
pca_at <- function(s) {
  if (is.null(qc) || is.null(qc$pcaData)) return(NULL)
  i <- match(s, qc$pcaData$name); if (is.na(i)) return(NULL)
  c(qc$pcaData$PC1[i], qc$pcaData$PC2[i])
}
xy_txt <- function(s) {
  v <- pca_at(s)
  if (is.null(v)) sprintf("%s(未在 PCA 坐标表中)", s) else sprintf("%s(%+.1f, %+.1f)", s, v[1], v[2])
}

## 3) 差异表达 (04): 上下调数、OE/RNAi 相关性、严格反向交集
o <- d$t_OE_CK; r <- d$t_RNAi_CK
union_deg <- union(union(d$deg_OE, d$deg_RNAi), d$deg_OR)
cor_oe_ri <- cor(o$log2FoldChange, r$log2FoldChange[match(o$gene_id, r$gene_id)],
                 use = "complete.obs")
deg_up <- function(set, tab) sum(set %in% tab$gene_id[!is.na(tab$log2FoldChange) & tab$log2FoldChange >  0])
deg_dn <- function(set, tab) sum(set %in% tab$gene_id[!is.na(tab$log2FoldChange) & tab$log2FoldChange <  0])
oe_up_n <- deg_up(d$deg_OE, o);   oe_dn_n <- deg_dn(d$deg_OE, o)
ri_up_n <- deg_up(d$deg_RNAi, r); ri_dn_n <- deg_dn(d$deg_RNAi, r)
sig_cell <- function(tab, dir) tab$gene_id[!is.na(tab$padj) & tab$padj < thr_p &
                                           !is.na(tab$log2FoldChange) &
                                           dir * tab$log2FoldChange >= thr_l]
strict_act <- length(intersect(sig_cell(o, 1),  sig_cell(r, -1)))
strict_rep <- length(intersect(sig_cell(o, -1), sig_cell(r, 1)))

## 4) 最强 DEG: 三次对比并集中, 在 OE 表里 padj 最小者
mm <- match(union_deg, o$gene_id)
top_gene <- if (length(union_deg)) union_deg[which.min(o$padj[mm])] else NA_character_
have_top <- !is.na(top_gene) && length(top_gene) == 1
if (have_top) {
  top_oe <- o[match(top_gene, o$gene_id), ]
  top_ri <- r[match(top_gene, r$gene_id), ]
  gcnt <- function(g) f$count_mat[top_gene, sample_cols[cond == g]]
  cnt_txt <- function(g) sprintf("%s (均值 %.0f)", paste(gcnt(g), collapse = ", "), mean(gcnt(g)))
  ds <- if (is.na(top_oe$gene_description)) "" else top_oe$gene_description
  ds <- trimws(strsplit(ds, "&&", fixed = TRUE)[[1]][1])
  if (nchar(ds) > 70) ds <- paste0(substr(ds, 1, 70), "…")
  same_dir_top <- !is.na(top_oe$log2FoldChange) && !is.na(top_ri$log2FoldChange) &&
                  sign(top_oe$log2FoldChange) == sign(top_ri$log2FoldChange)
}

## 5) 富集 (07): 逐表的最小 p / 最小 padj / 显著数; GO 与 KEGG 注释规模
en_rds <- file.path(rdsdir, "07_enrichment.rds")
en <- if (file.exists(en_rds)) readRDS(en_rds) else NULL
einfo <- function(key) {
  e <- if (is.null(en)) NULL else en$enr[[key]]
  if (is.null(e))     return(list(n = NA_integer_, nsig = NA_integer_, minp = NA_real_,
                                  minpadj = NA_real_, feat = NA_character_))
  if (nrow(e) == 0)   return(list(n = 0L, nsig = 0L, minp = NA_real_,
                                  minpadj = NA_real_, feat = NA_character_))
  i <- which.min(e$pvalue)
  list(n = nrow(e), nsig = sum(e$padj < 0.05, na.rm = TRUE),
       minp = min(e$pvalue, na.rm = TRUE), minpadj = min(e$padj, na.rm = TRUE),
       feat = as.character(e$feature[i]))
}
fmt_sci  <- function(x, d = 1) if (is.na(x)) "—" else formatC(x, format = "e", digits = d)
fmt_padj <- function(x)          if (is.na(x)) "—" else formatC(x, format = "g", digits = 3)
enr_row  <- function(label, key) {
  i <- einfo(key)
  sprintf("| %s | %s | %s | %s |", label, fmt_sci(i$minp), fmt_padj(i$minpadj),
          if (is.na(i$nsig)) "未产出" else if (i$nsig > 0) "是" else "否")
}
enr_feat_txt <- function(key) {
  i <- einfo(key)
  if (is.na(i$n)) sprintf("未产出 (%s 表缺失)", key)
  else if (i$n == 0) sprintf("空集 (测试集为空, 无检验)")
  else sprintf("%s (最小 raw p %s, BH padj %s)", i$feat, fmt_sci(i$minp), fmt_padj(i$minpadj))
}
go_stats   <- if (is.null(en)) NULL else en$go_stats
kegg_stats <- if (is.null(en)) NULL else en$kegg_stats
go_ann <- if (file.exists(file.path(rdsdir, "06_go.rds"))) readRDS(file.path(rdsdir, "06_go.rds")) else NULL
go_ok  <- !is.null(go_stats) ||
          (!is.null(go_ann) && !is.null(go_ann$go_term2gene) && nrow(go_ann$go_term2gene) > 0)
n_go_genes <- if (!is.null(go_stats)) go_stats$genes else
              if (!is.null(go_ann)) length(unique(go_ann$go_term2gene$gene_id)) else NA_integer_
n_go_terms <- if (!is.null(go_stats)) go_stats$terms else
              if (!is.null(go_ann)) length(unique(go_ann$go_term2gene$go_id)) else NA_integer_
go_per_gene <- if (!is.null(go_stats)) round(go_stats$gene_term / go_stats$genes, 1) else NA_real_
go_density_txt <- if (is.na(go_per_gene)) "注释密度未知" else sprintf("平均 %.0f GO/基因", go_per_gene)
go_src <- if (!is.null(go_stats)) go_stats$source else if (!is.null(go_ann)) go_ann$source else NA_character_
go_src_txt <- if (!go_ok) {
  "**GO 注释未产出** (06 步骤失败/未运行)"
} else if (identical(go_src, "local")) {
  sprintf("**目标物种原生(eggNOG)本地注释**:%s 基因 / %s 个 GO",
          format(n_go_genes, big.mark = ","), format(n_go_terms, big.mark = ","))
} else {
  sprintf("UniProt 同源转移注释:%s 基因 / %s 个 GO",
          format(n_go_genes, big.mark = ","), format(n_go_terms, big.mark = ","))
}
heat_top <- as.integer(Sys.getenv("HEATMAP_TOP", "100"))

## 9) 富集背景口径 (有资格判 DEG 的基因) 与敏感性对照
bgs  <- if (is.null(en)) NULL else en$backgrounds
sens <- if (is.null(en)) NULL else en$sens
bg_txt <- if (!is.null(bgs)) {
  sprintf("背景为**有资格被判为 DEG 的基因**(逐分析取 padj 非 NA 者:OE_vs_CK %d、RNAi_vs_CK %d、OE_vs_RNAi %d、并集 %d、靶标 %d;全部过滤基因为 %d)",
          length(bgs$OE_vs_CK), length(bgs$RNAi_vs_CK), length(bgs$OE_vs_RNAi),
          length(bgs$union), length(bgs$target), length(bgs$all_filtered))
} else {
  "背景定义未记录(07 未产出 backgrounds 字段)"
}
n_sens_diff <- if (!is.null(sens) && nrow(sens)) sum(sens$n_sig_tested_bg != sens$n_sig_allfiltered_bg) else NA_integer_
sens_txt <- if (is.na(n_sens_diff)) {
  "**背景口径敏感性**:未做(缺少对照结果)。"
} else {
  sprintf("**背景口径敏感性**:同一批检验改用“全部过滤基因”背景复算,padj<0.05 项数在 %d/%d 张表上完全一致%s(逐表对照见 `tables/enrich_background_sensitivity.csv`)。",
          nrow(sens) - n_sens_diff, nrow(sens),
          if (n_sens_diff > 0) sprintf(",%d 张表不同", n_sens_diff) else "")
}

## 6) 注释检索: 目标基因能否按名定位
anno <- f$anno
target_symbol <- Sys.getenv("TARGET_GENE_SYMBOL", "TF")
n_tf_named <- sum(grepl(target_symbol, anno$gene_name, ignore.case = TRUE), na.rm = TRUE) +
              sum(grepl(target_symbol, anno$gene_description, ignore.case = TRUE), na.rm = TRUE)

## 7) 火山图"高精度小效应"示例: OE 表中 |LFC| < 阈值却 padj 最小的基因
tiny_i <- which(!is.na(o$padj) & abs(o$log2FoldChange) < thr_l)
tiny_txt <- if (length(tiny_i) > 0) {
  k <- tiny_i[which.min(o$padj[tiny_i])]
  sprintf("`%s`(log2FC≈**%+.3f**,padj≈%s)", o$gene_id[k], o$log2FoldChange[k],
          fmt_sci(o$padj[k], 2))
} else "(本次未出现该情况)"

## 8) PCA 组内离散 vs 组间距离 (用于判断处理效应是否盖过样本异质性)
pca_within_txt <- "(PCA 坐标缺失)"; pca_between_txt <- ""
pca_verdict <- "样本离散度无法评估(PCA 坐标缺失)。"
if (!is.null(qc) && !is.null(qc$pcaData)) {
  pd <- qc$pcaData
  grps <- c("CK", "RNAi", "OE")
  within_max <- sapply(grps, function(g) {
    s <- pd[pd$condition == g, c("PC1", "PC2"), drop = FALSE]
    if (nrow(s) < 2) NA_real_ else max(dist(s))
  })
  cen <- t(sapply(grps, function(g) colMeans(pd[pd$condition == g, c("PC1", "PC2"), drop = FALSE])))
  btw <- as.matrix(dist(cen)); bv <- btw[upper.tri(btw)]
  pca_within_txt <- paste(sprintf("%s=%.1f", grps, within_max), collapse = ";")
  pca_between_txt <- sprintf("组间质心距离:CK–RNAi=%.1f、CK–OE=%.1f、RNAi–OE=%.1f。",
                             btw["CK","RNAi"], btw["CK","OE"], btw["RNAi","OE"])
  pca_verdict <- if (max(within_max, na.rm = TRUE) >= min(bv)) {
    "⚠ 组内离散已达到组间距离量级,提示样本异质性/批次效应不可忽略:严格靶标与富集结论需谨慎(见第 8 节)。"
  } else {
    "组内离散小于组间距离,处理效应方向可用。"
  }
}

L <- c(
'# 目标转录因子 RNA-seq 分析报告','',
'**物种**:目标物种（参考基因组）',
'**数据**:目标 TF 过表达 / 敲低 / 对照转录组',
'**流程**:`Pipeline/01…10`(分步运行, DESeq2 全流程);本报告数字全部由产物文件现算',
sprintf('**生成时间**:%s ・ 所有数字可在 `results/` 对应文件中核对', Sys.time()),'',
'---','',
'## 0. 摘要(先看这里)','',
'| 关键指标 | 结果 | 出处 |','|---|---|---|',
sprintf('| 原始 / 过滤后基因 | %s → **%s** | `tables/filter_stats.csv` |',
        format(n_raw, big.mark=","), format(n_gene, big.mark=",")),
sprintf('| DEG OE vs CK / RNAi vs CK / OE vs RNAi | **%s / %s / %s** | `summary.csv` |', n_oe, n_ri, n_or),
sprintf('| 最强且最稳 DEG | `%s`%s | `tables/DEG_*.csv` |',
        if (have_top) top_gene else "(无 DEG)",
        if (have_top) sprintf('(%s),%s', top_oe$Family,
                              if (same_dir_top) "OE 与 RNAi 中**同向**变化" else "OE 与 RNAi 中**反向**变化") else ""),
sprintf('| 严格反向靶标(激活/抑制) | **%d / %d** | `summary.csv` |', n_act, n_rep),
sprintf('| 探索性宽松靶标(激活/抑制) | **%d / %d** | `tables/targets_*_exploratory.csv` |', n_ract, n_rrep),
sprintf('| GO / PFAM / TF 富集显著数 (BH<0.05) | GO=%s, KEGG=%s, PFAM=%s, TF=%s | `tables/enrich_*.csv` |',
        ifelse(is.na(einfo("GO_DEG_union")$nsig), "未产出", einfo("GO_DEG_union")$nsig),
        ifelse(is.na(einfo("KEGG_DEG_union")$nsig), "未产出", einfo("KEGG_DEG_union")$nsig),
        ifelse(is.na(einfo("PFAM_DEG_union")$nsig), "未产出", einfo("PFAM_DEG_union")$nsig),
        ifelse(is.na(einfo("TFfamily_DEG_union")$nsig), "未产出", einfo("TFfamily_DEG_union")$nsig)),
sprintf('| 高分候选位点 (OE_only, score 最高) | `%s` | `tables/candidate_TF_locus.csv` |',
        if (nrow(tg$cand_TF) > 0) tg$cand_TF$gene_id[1] else "(无)"),
sprintf('| 目标基因本体 (%s) | 注释按名检索 "%s" 命中 **%d** 行 → %s | `tables/candidate_TF_locus.csv` |',
        target_symbol, target_symbol, n_tf_named,
        if (n_tf_named == 0) "**无法按名定位**,身份需 BLAST 确认" else "存在同名/同源注释,仍建议序列确认"),'',
sprintf('**一句话结论**:过表达与敲低两组扰动的 log2FC 相关系数为 **%.3f**,严格双向靶标 **%d** 个;%s PFAM/TF 家族富集显著 %d 项。%s',
        cor_oe_ri, n_act + n_rep,
        if (!go_ok) "⚠ 本次**未产出 GO 富集**(06 步骤未给出可用注释, 故无从判断 GO 是否显著)。"
        else if (einfo("GO_DEG_union")$nsig == 0)
          sprintf("GO 富集 0 项显著(%s, 密集注释稀释信号),", go_density_txt)
        else sprintf("GO 富集 %d 项显著,", einfo("GO_DEG_union")$nsig),
        einfo("PFAM_DEG_union")$nsig + einfo("TFfamily_DEG_union")$nsig,
        "稳健结论以 QC / DE / 靶标方向等 count 层面观察为准;探索性候选与候选位点需 qRT-PCR / BLAST 验证。"),'',
'---','',
'## 1. 实验设计与样本结构','',
sprintf('- 分组:**OE(%d)、RNAi(%d)、CK(%d)**;设计 `%s`(CK 为参照)。',
        sum(cond == "OE"), sum(cond == "RNAi"), sum(cond == "CK"),
        if (!is.null(d$params$design)) d$params$design else "~ condition"),
sprintf('- 低表达过滤:count >= %s 且至少 %s 个样本 (见 `tables/preprocessing_methods.txt`)。',
        f_min_count, f_min_samples),
sprintf('- %s **CK 仅 %d 个生物学重复**:DESeq2 离散度估计不稳、假阳性偏高,流程已就此**明确告警**。此为实验设计局限,非代码问题。',
        if (sum(cond == "CK") < 3) "⚠️" else "✅", sum(cond == "CK")),'',
'---','',
'## 2. 质控:主成分分析(PCA,步骤 03,置于 DE 前)','',
'**图**:`figures/PCA.svg` ・ 坐标:`tables/PCA_coords.csv` ・ 方法:`tables/viz_methods.txt`',
sprintf('方法:对 DESeq2 `vst`(方差稳定转换,blind=TRUE)后的表达矩阵,取方差最高的基因做 PCA。'),'',
'| 主成分 | 方差贡献 |','|---|---|',
sprintf('| PC1 | **%s** |', pct_txt(1)),
sprintf('| PC2 | **%s** |', pct_txt(2)),
sprintf('| PC3 | %s |', pct_txt(3)),'',
'**读图要点(均可由 PCA_coords.csv 复算)**:',
sprintf('- 两个 CK:%s、%s。', xy_txt("CK1"), xy_txt("CK2")),
sprintf('- OE 组内:%s、%s、%s;RNAi 组内:%s、%s、%s。',
        xy_txt("OE3"), xy_txt("OE1"), xy_txt("OE2"),
        xy_txt("KD1"), xy_txt("KD2"), xy_txt("KD3")),
sprintf('- 组内最大两两距离(PC1–PC2):%s。%s', pca_within_txt, pca_between_txt),
sprintf('- 结论:%s', pca_verdict),'',
'---','',
'## 3. 差异表达分析(DESeq2,步骤 04)','',
sprintf('**方法/文献**:DESeq2(Love et al., 2014);log2FC 用 `lfcShrink`(type=%s;apeglm: Zhu et al., 2019)收缩;padj 为 Benjamini–Hochberg FDR(1995);独立过滤(Bourgon et al., 2010)。阈值 **padj<%g 且 |log2FC|≥%g**。详见 `tables/DEG_methods.txt`。', shrink_type, thr_p, thr_l),'',
'| 对比 | DEG 数 | 上调 / 下调 | 图 |','|---|---|---|---|',
sprintf('| OE vs CK | %d | %d / %d | `figures/Volcano_OE_vs_CK.svg` |', n_oe, oe_up_n, oe_dn_n),
sprintf('| RNAi vs CK | %d | %d / %d | `figures/Volcano_RNAi_vs_CK.svg` |', n_ri, ri_up_n, ri_dn_n),
sprintf('| OE vs RNAi | %d | — | `figures/Volcano_OE_vs_RNAi.svg` |', n_or),
'| 三组重叠 | — | — | `figures/Venn_DEG.svg` |','',
sprintf('**⚠️ 解读火山图 / DEG 表的重要提醒**:`DEG_*.csv` 按 padj 排序的"头部"里混有**极低 |log2FC| 但极高精度**的基因(如 %s)——它们不满足 |log2FC|≥%g,**不计入 DEG**。这正是必须用 "padj **且** |LFC|" 双阈值的原因。',
        tiny_txt, thr_l),'',
if (have_top) sprintf('**最强 DEG(并集中 OE 表 padj 最小者)**:`%s`(Family=%s;注释首段:%s),baseMean≈%.0f:',
                      top_gene, top_oe$Family, ds, top_oe$baseMean) else '**(无 DEG)**',
if (have_top) c('','| 分组 | 重复原始计数 | DESeq2 vs CK |','|---|---|---|',
  sprintf('| RNAi | %s | log2FC **%+.2f**,padj %s |', cnt_txt("RNAi"), top_ri$log2FoldChange, fmt_sci(top_ri$padj, 2)),
  sprintf('| OE | %s | log2FC **%+.2f**,padj %s |', cnt_txt("OE"), top_oe$log2FoldChange, fmt_sci(top_oe$padj, 2)),
  sprintf('| CK | %s | 参照 |', cnt_txt("CK")), '') else NULL,
if (have_top) sprintf('它在 OE 与 RNAi 中%s——%s', if (same_dir_top) "**同向**变化" else "**反向**变化",
  if (same_dir_top) "更像对两种扰动/胁迫普遍响应的基因,**不是** TF 的干净反向靶标。"
  else "符合反向调控模型,可作为反向靶标的优先验证对象。") else NULL,'',
'---','',
'## 4. TF 靶标鉴定(核心,步骤 05)','',
'**策略**:利用过表达/敲低的**反向变化**筛选靶标——激活=OE 显著上调**且** RNAi 显著下调;抑制=反向。这是 TF 靶标鉴定的**标准保守遗传策略**(gain/loss-of-function 反向验证)。',
'**图**:`figures/Scatter_target_direction.svg`。','',
'| 靶标集 | 激活 | 抑制 |','|---|---|---|',
sprintf('| 严格(双侧显著反向) | **%d** | **%d** |', n_act, n_rep),
sprintf('| 探索性(一侧显著+另一侧反向趋势) | **%d** | **%d** |', n_ract, n_rrep),'',
sprintf('**为何严格靶标是 %d**:OE 上调 %d / 下调 %d;RNAi 上调 %d / 下调 %d;`OE↑ ∩ RNAi↓ = %d`、`OE↓ ∩ RNAi↑ = %d`;OE 与 RNAi 的 log2FC 相关系数 %.3f。',
        n_act, oe_up_n, oe_dn_n, ri_up_n, ri_dn_n, strict_act, strict_rep, cor_oe_ri),'',
'**候选位点(命名与证据澄清)**:

- score = log2FC(OE/CK) − log2FC(RNAi/CK) 是**项目自定义启发式,非公认统计量**。它按"哪一侧幅度大"排序,会把两类不同基因排在一起:**仅 OE 侧显著变化**(过表达特异/转基因插入位点候选)与**仅 RNAi 侧显著变化**(敲低特异候选靶标)。
- 因此 `tables/candidate_reciprocal_score.csv` 给出全部 score top-N 及其**双侧 log2FC/padj、baseMean 和 evidence 分级**(OE_only / RNAi_only / both_signif_opposite / both_signif_same_direction / neither_signif),可逐条核验;`tables/candidate_TF_locus.csv` 只保留 **OE_only** 子集,才与"TF 过表达位点候选"的命名相符。
- ⚠ 该 score 本身不含证据门槛(不看显著性、不看表达量),排序靠前不等于可信;位点身份仍需 BLAST 确认。此前的混合 top-N 表已停用。','',
sprintf('- 本次 OE_only 候选 %d 个(score 最高:%s);若同时需要"敲低特异"候选,见 `candidate_reciprocal_score.csv` 中 evidence=RNAi_only 的行。',
        nrow(tg$cand_TF), if (nrow(tg$cand_TF) > 0) tg$cand_TF$gene_id[1] else "(无)"),'',
sprintf('- **目标基因身份**:注释表 gene_name / gene_description 中检索 "TF" 命中 %d 行 → %s。',
        n_tf_named,
        if (n_tf_named == 0)
          "**无法按名定位**(gene_name 为基因模型 ID),目标基因身份必须用序列(BLAST)确认"
        else "存在同名/同源注释,仍建议用序列确认"),'',
'---','',
'## 5. 功能富集分析(步骤 06 GO 注释 + 07 富集)','',
sprintf('**方法**:GO 注释为 %s(KEGG 通路用目标物种原生 KO%s);PFAM/TF 家族用超几何检验,特异性用 2×2 Fisher,均 BH 校正(1995)。%s —— 不能判 DEG 的基因(padj=NA,多为极低表达)不进入背景,否则会破坏超几何检验的抽样模型;方法与文献见 `tables/enrichment_methods.txt`。**图**:`figures/Bubble_GO_*.svg`、`figures/Bubble_KEGG_*.svg`;x=GeneRatio、点大小=Count、颜色=-log10(padj)。',
  go_src_txt,
  if (!is.null(kegg_stats)) sprintf("(%s 基因 / %s KO)", format(kegg_stats$genes, big.mark=","),
                                    format(kegg_stats$kos, big.mark=",")) else "(KEGG 数据缺失, 未做 KEGG 富集)",
  bg_txt),'',
sprintf('**结果**:%sPFAM/TF 家族富集显著 %d 项。',
  if (!go_ok) "⚠ 本次**未产出 GO 富集**(06 步骤未给出可用注释;这不等于“无富集”)。"
  else if (all(c(einfo("GO_DEG_union")$nsig, einfo("GO_OE")$nsig, einfo("GO_RNAi")$nsig,
                  einfo("GO_OE_RNAi")$nsig) == 0)) "GO/KEGG 富集:全部基因集 0 项 BH<0.05 显著。"
  else sprintf("GO 显著性见下表(padj<0.05:%d/%d/%d/%d)。", einfo("GO_DEG_union")$nsig,
               einfo("GO_OE")$nsig, einfo("GO_RNAi")$nsig, einfo("GO_OE_RNAi")$nsig),
  einfo("PFAM_DEG_union")$nsig + einfo("TFfamily_DEG_union")$nsig + einfo("PFAM_spec")$nsig +
  einfo("TFfamily_spec")$nsig),'',
'| 基因集 | 最高 raw p | BH padj | 是否显著 (padj<0.05) |','|---|---|---|---|',
enr_row('GO(DEG 并集)',    'GO_DEG_union'),
enr_row('GO(OE vs CK)',    'GO_OE'),
enr_row('GO(RNAi vs CK)',  'GO_RNAi'),
enr_row('GO(OE vs RNAi)',  'GO_OE_RNAi'),
enr_row('GO(激活靶标)',    'GO_activated'),
enr_row('GO(抑制靶标)',    'GO_repressed'),
enr_row('KEGG(DEG 并集)',  'KEGG_DEG_union'),
enr_row('KEGG(OE vs CK)',  'KEGG_OE'),
enr_row('KEGG(RNAi vs CK)','KEGG_RNAi'),
enr_row('KEGG(OE vs RNAi)','KEGG_OE_RNAi'),
enr_row('PFAM(DEG 并集)',  'PFAM_DEG_union'),
enr_row('TF 家族(DEG 并集)','TFfamily_DEG_union'),
'','**各集最小 p 的特征**:',
sprintf('- GO(DEG 并集):%s', enr_feat_txt('GO_DEG_union')),
sprintf('- GO(OE vs CK):%s', enr_feat_txt('GO_OE')),
sprintf('- GO(RNAi vs CK):%s', enr_feat_txt('GO_RNAi')),
sprintf('- GO(OE vs RNAi):%s', enr_feat_txt('GO_OE_RNAi')),
sprintf('- KEGG(DEG 并集):%s', enr_feat_txt('KEGG_DEG_union')),
sprintf('- PFAM(DEG 并集):%s', enr_feat_txt('PFAM_DEG_union')),
sprintf('- TF 家族(DEG 并集):%s', enr_feat_txt('TFfamily_DEG_union')),
'', sens_txt, '',
sprintf('**为何可能 0 显著**:(1) 注释密度高时(本次%s)背景计数偏高、超几何被稀释;(2) 严格靶标集为空 -> 对应富集无测试集;(3) 背景大、测试集小 -> 检验效能有限。', go_density_txt),
'**注释来源敏感性**:不同注释来源的覆盖密度不同，GO 富集结论也可能随之变化。流程会记录注释来源、版本与匹配率，并并排输出两种背景口径的敏感性对照，便于判断结论是否稳健。','',
'**结论与建议**:',
'- 稳健结论以 count 层面的 DE / 靶标方向 / QC 观察为准;相对可靠的注释证据是 PFAM/TF 家族(来自序列结构域)。',
'- 若需 GO 富集,建议对密集注释**去冗余**(只留较特异 GO,或用 topGO/clusterProfiler 的 elim/weight 按 GO 树层级校正)后再试。',
'- 详见 `README.md`（本地 GO 注释工作流与密度权衡）。','',
'---','',
'## 6. 热图(Top DEG 表达模式,步骤 08)','',
sprintf('**图**:`figures/Heatmap_topDEG.svg` ・ 取三次对比 DEG **并集**按 |log2FC(OE/CK)| 降序前 **%d** 个(HEATMAP_TOP),log2(count+1) 后按基因行 z-score 绘制。', heat_top),'',
'---','',
'## 7. 主要结论','',
sprintf('1. **样本结构**:PC1=%s、PC2=%s、PC3=%s;%s',
        pct_txt(1), pct_txt(2), pct_txt(3),
        '各组内离散与组间质心距离见第 2 节(可直接从 `tables/PCA_coords.csv` 复算)。'),
sprintf('2. **DEG 丰富**(OE %d / RNAi %d / OE vs RNAi %d)%s。', n_oe, n_ri, n_or,
        if (have_top) sprintf(',最强者为 `%s`(%s)且在两处理中%s', top_gene, top_oe$Family,
                              if (same_dir_top) "同向变化" else "反向变化") else ""),
sprintf('3. **严格反向靶标 %d 个**(OE/RNAi log2FC 相关 %.3f):%s',
        n_act + n_rep, cor_oe_ri,
        if (n_act + n_rep == 0) "该结果为真实结果(非 bug),提示 RNAi 敲低效应弱或组内离散大。" else "见 `tables/targets_activated.csv` / `targets_repressed.csv`。"),
sprintf('4. **富集**:%s PFAM/TF 家族显著 %d 项;%s',
        if (!go_ok) "GO 未产出(见第 5 节);" else if (einfo("GO_DEG_union")$nsig == 0) "GO 0 项显著;" else
          sprintf("GO %d 项显著;", einfo("GO_DEG_union")$nsig),
        einfo("PFAM_DEG_union")$nsig + einfo("TFfamily_DEG_union")$nsig,
        if (!go_ok) "不要将 GO 缺失解读为无富集。" else "富集结论受注释密度与效能限制。"),
sprintf('5. **高分候选位点**:`%s`(需 BLAST 确认身份;score 为项目自定义启发式,非公认统计量)。',
        if (nrow(tg$cand_TF) > 0) tg$cand_TF$gene_id[1] else "无"),'',
'---','',
'## 8. 局限与注意事项','',
sprintf('- **CK n=%d**:参照组离散度估计不稳,假阳性风险偏高(已告警)。', sum(cond == "CK")),
'- **批次效应未校正**:当前设计见 `tables/DEG_methods.txt`;如测序/建库存在批次,可设 `DESEQ_BATCH="..."` 切成 `~ batch + condition` 后重跑 04(注意批次与分组完全混杂时设计矩阵不可解)。',
sprintf('- **GO 注释**:%s。注释来源与密度直接决定富集结论,换来源前不要比较 padj。', go_src_txt),
'- **富集效能有限**:背景大、测试集小;需去冗余算法或更多重复才能提升。',
sprintf('- **候选 score 无证据门槛**:排序不能替代显著性判断,`candidate_reciprocal_score.csv` 中已给出双侧 padj/baseMean 供核验。'),
'','---','',
'## 9. 建议的下一步','',
sprintf('1. **验证候选位点身份**:对 `%s` 等候选做 BLAST;并确认 RNAi 是否真的降低了目标转录本。',
        if (nrow(tg$cand_TF) > 0) tg$cand_TF$gene_id[1] else "候选位点"),
'2. **批次校正**:设 `DESEQ_BATCH="b1,b1,..."`(与样本列同序)后重跑 04→07,与当前结果对比稳健性。',
'3. **增加重复**:尤其 CK 与 RNAi,提升离散度估计与富集效能。',
sprintf('4. **qRT-PCR 验证**:探索性宽松靶标(激活 %d / 抑制 %d)中的高置信候选 + 效应最强的 DEG。', n_ract, n_rrep),
'5. **如需 GO 富集**:对密集注释去冗余(topGO elim / clusterProfiler)后重跑 06→07。','',
'---','',
'## 10. 可复现性与文献','',
'**分步运行**(本流程不提供一键包装脚本):','','```bash',
'for s in 01_read_data 02_filter 03_qc 04_deseq 05_targets 06_go 07_enrichment 08_visualization 09_summary 10_report; do',
'  Rscript Pipeline/$s.R || break',
'done','```','',
sprintf('- 过滤/QC/DE 阈值经环境变量覆盖(本次: `FILTER_MIN_COUNT=%s FILTER_MIN_SAMPLES=%s DEG_PADJ=%g DEG_LFC=%g DESEQ_SHRINK=%s`)',
        f_min_count, f_min_samples, thr_p, thr_l, Sys.getenv("DESEQ_SHRINK", "apeglm")),
sprintf('- **GO 来源**:%s;设 `GO_ANNOT_FILE=...` 可指定其它本地文件(`run_eggnog.sh` + `build_gene_go.R` 生成)。', go_src_txt),
'- **方法文档**:`tables/preprocessing_methods.txt`、`tables/DEG_methods.txt`、`tables/targets_methods.txt`、`tables/enrichment_methods.txt`(含背景口径与文献)、`tables/viz_methods.txt`。',
'- **稳健性对照**:`tables/enrich_background_sensitivity.csv`(两种背景口径逐表对照)。',
'- **来源留痕**:`tables/sessionInfo.txt`(R 与包版本)、`tables/GO_source_provenance.txt`、`tables/KEGG_source_provenance.txt`(注释/数据库文件路径 + 大小 + md5 + 时间)。',
'- **主文献清单(含 DOI/PMID)**:`tables/references.md`;**目标物种 GO 工作流**:`README.md`。',
'- 关键文献:DESeq2(Love 2014)、apeglm(Zhu 2019)、BH(1995)、limma(Ritchie 2015)、Eisen 1998、eggNOG-mapper(Huerta-Cepas 2017)、PlantTFDB(Jin 2017)。','',
'---','',
'*本报告由 `Pipeline/10_report.R` 生成;所有数字在生成时从产物文件现算。凡自定义/非标准步骤(target score、宽松靶标、OE_only 候选分级)均已标注,不与公认统计量混淆。*'
)
L <- L[!vapply(L, is.null, logical(1))]
writeLines(L, file.path(outdir, "ANALYSIS_REPORT.md"))

## 运行环境留痕 (论文可复现性: 记录 R 与关键包版本)
##   TZ 未设置时 sessionInfo() 会去调用 timedatectl, 在无 systemd 的容器里会抛
##   "Failed to connect to bus" 到 stderr(只是噪音); 先固定 TZ 让日志干净可复现。
tz_env <- Sys.getenv("TZ")
if (!nzchar(tz_env)) Sys.setenv(TZ = "Asia/Shanghai")
writeLines(capture.output(suppressWarnings(suppressMessages(sessionInfo()))),
           file.path(outdir, "tables", "sessionInfo.txt"))
if (!nzchar(tz_env)) Sys.unsetenv("TZ")

cat(sprintf("[%s] 分析报告已生成: %s (%d 行; 数字均为现算)\n",
            Sys.time(), file.path(outdir, "ANALYSIS_REPORT.md"), length(L)))
