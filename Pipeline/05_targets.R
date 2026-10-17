# ============================================================
# 05_targets.R  —  OE/RNAi 反向一致性 -> TF 候选靶基因
#   激活靶标: OE显著上调 且 RNAi显著下调;  抑制靶标: OE显著下调 且 RNAi显著上调
#   (转录组靶标常用保守策略: 双侧反向一致变化 -> 直接/间接靶标)
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

d <- readRDS(file.path(rdsdir, "04_deseq.rds"))
t_OE_CK   <- d$t_OE_CK
t_RNAi_CK <- d$t_RNAi_CK
anno      <- readRDS(file.path(rdsdir, "02_filtered.rds"))$anno

## 3) 阈值复用 04 保存的参数 (单一事实来源, 改 04 阈值时此处自动同步)
if (is.null(d$params)) { padj_thresh <- 0.05; lfc_thresh <- 1
} else { padj_thresh <- d$params$padj_thresh; lfc_thresh <- d$params$lfc_thresh }
top_n <- as.integer(Sys.getenv("TARGET_TOP_N", "10"))

## 2) 兼容校验: 两张差异表基因集须一致 (顺序可不同, 下面按基因名对齐)
if (!setequal(t_OE_CK$gene_id, t_RNAi_CK$gene_id))
  stop("t_OE_CK 与 t_RNAi_CK 基因集合不一致, 请检查 03 输出。")

oe_lfc  <- setNames(t_OE_CK$log2FoldChange,   t_OE_CK$gene_id)
oe_padj <- setNames(t_OE_CK$padj,             t_OE_CK$gene_id)
ri_lfc  <- setNames(t_RNAi_CK$log2FoldChange, t_RNAi_CK$gene_id)
ri_padj <- setNames(t_RNAi_CK$padj,           t_RNAi_CK$gene_id)

sig_up   <- function(lfc, padj) !is.na(lfc) & !is.na(padj) & padj < padj_thresh & lfc >=  lfc_thresh
sig_down <- function(lfc, padj) !is.na(lfc) & !is.na(padj) & padj < padj_thresh & lfc <= -lfc_thresh

## 1) 严格保守筛选: 激活=OE上调 且 RNAi下调; 抑制=OE下调 且 RNAi上调
##    (双侧均显著以控制假阳性; 该规则为方法学既定设计, 不放松)
## 按基因名(而非位置)提取显著基因集, 避免两张表顺序不同导致的错位
oe_up   <- names(oe_lfc)[sig_up(oe_lfc, oe_padj)]
ri_down <- names(ri_lfc)[sig_down(ri_lfc, ri_padj)]
oe_down <- names(oe_lfc)[sig_down(oe_lfc, oe_padj)]
ri_up   <- names(ri_lfc)[sig_up(ri_lfc, ri_padj)]
act <- intersect(oe_up, ri_down)
rep <- intersect(oe_down, ri_up)

## 6) 校验: 靶标基因须能在 anno 中找到 (避免 match 输出 NA 行)
if (!all(act %in% anno$gene_id)) stop("部分激活靶标 gene_id 不在 anno 中")
if (!all(rep %in% anno$gene_id)) stop("部分抑制靶标 gene_id 不在 anno 中")
act_anno <- anno[match(act, anno$gene_id), ]; rownames(act_anno) <- NULL
rep_anno <- anno[match(rep, anno$gene_id), ]; rownames(rep_anno) <- NULL
write.csv(act_anno, file.path(outdir, "tables", "targets_activated.csv"),   row.names = FALSE)
write.csv(rep_anno, file.path(outdir, "tables", "targets_repressed.csv"),  row.names = FALSE)

## 1) 探索性(宽松)候选: 至少一侧显著 且 另一侧呈一致反向趋势 (不要求双侧均显著)
##    明确标记为 exploratory —— 用于弥补严格双侧阈值可能漏掉的趋势性候选, 不作主结论
ri_down_trend <- names(ri_lfc)[ri_lfc < 0]   # RNAi 呈下调趋势(不限显著性)
ri_up_trend   <- names(ri_lfc)[ri_lfc > 0]
oe_down_trend <- names(oe_lfc)[oe_lfc < 0]
oe_up_trend   <- names(oe_lfc)[oe_lfc > 0]
relaxed_act <- union(oe_up[oe_up %in% ri_down_trend],
                     ri_down[ri_down %in% oe_up_trend])
relaxed_rep <- union(oe_down[oe_down %in% ri_up_trend],
                     ri_up[ri_up %in% oe_down_trend])
if (length(relaxed_act) && all(relaxed_act %in% anno$gene_id))
  write.csv(anno[match(relaxed_act, anno$gene_id), ],
            file.path(outdir, "tables", "targets_activated_exploratory.csv"), row.names = FALSE)
if (length(relaxed_rep) && all(relaxed_rep %in% anno$gene_id))
  write.csv(anno[match(relaxed_rep, anno$gene_id), ],
            file.path(outdir, "tables", "targets_repressed_exploratory.csv"), row.names = FALSE)

## 候选位点打分 (项目自定义启发式, 非公认统计量; 详见 targets_methods.txt)
## 4) score = log2FC(OE/CK) - log2FC(RNAi/CK): 真实靶标两侧反向, 差值=两效应绝对值之和。
##    ⚠ 该 score 会把"两类完全不同的基因"排到前面, 故必须分表输出并附证据列, 不能只输出
##      一个混合 top-N 就当作 "TF 位点候选":
##        * OE_only    : OE 显著变化而 RNAi 无显著变化 -> 过表达特异(转基因/插入位点候选)
##        * RNAi_only  : RNAi 显著变化而 OE 无显著变化 -> 敲低特异(候选靶标, 属另一类结论)
##        * 其余       : both_signif_* / neither_signif
##    输出:
##      candidate_reciprocal_score.csv = 全部基因按 score 降序 top-N (含证据列, 透明可核验)
##      candidate_TF_locus.csv      = 其中 OE_only 子集 top-N (才与"TF 过表达位点"命名相符)
ri_lfc_all  <- t_RNAi_CK$log2FoldChange[match(t_OE_CK$gene_id, t_RNAi_CK$gene_id)]
ri_padj_all <- t_RNAi_CK$padj[match(t_OE_CK$gene_id, t_RNAi_CK$gene_id)]
keep_cand   <- !is.na(t_OE_CK$log2FoldChange) & !is.na(ri_lfc_all)
cand <- t_OE_CK[keep_cand, ]
cand$log2FC_RNAi <- ri_lfc_all[keep_cand]
cand$padj_RNAi   <- ri_padj_all[keep_cand]
cand$score       <- cand$log2FoldChange - cand$log2FC_RNAi

## 证据分级: 至少一侧 padj<阈值 才算"显著", 双侧显著且方向相反即严格反向靶标
oe_sig <- !is.na(cand$padj)      & cand$padj      < padj_thresh &
          !is.na(cand$log2FoldChange) & abs(cand$log2FoldChange) >= lfc_thresh
ri_sig <- !is.na(cand$padj_RNAi) & cand$padj_RNAi < padj_thresh &
          !is.na(cand$log2FC_RNAi)   & abs(cand$log2FC_RNAi)   >= lfc_thresh
ev <- rep("neither_signif", nrow(cand))
ev[ri_sig & !oe_sig] <- "RNAi_only"
ev[oe_sig & !ri_sig] <- "OE_only"
both <- oe_sig & ri_sig
same_dir <- both & sign(cand$log2FoldChange) == sign(cand$log2FC_RNAi)
ev[both & !same_dir] <- "both_signif_opposite"
ev[same_dir]         <- "both_signif_same_direction"
cand$evidence <- ev

cand <- cand[order(-cand$score), ]
cand_cols <- c("gene_id", "gene_name", "Family", "baseMean",
               "log2FoldChange", "padj", "log2FC_RNAi", "padj_RNAi",
               "score", "evidence")
cand_reciprocal <- head(cand, top_n)[, cand_cols]
write.csv(cand_reciprocal,
          file.path(outdir, "tables", "candidate_reciprocal_score.csv"), row.names = FALSE)

cand_TF <- head(cand[cand$evidence == "OE_only", , drop = FALSE], top_n)[, cand_cols]
if (nrow(cand_TF) > 0) {
  write.csv(cand_TF,
            file.path(outdir, "tables", "candidate_TF_locus.csv"), row.names = FALSE)
} else {
  warning("无 OE_only 候选 (OE 显著变化且 RNAi 无显著变化的基因), 未写出 candidate_TF_locus.csv")
}

saveRDS(list(act = act, rep = rep, cand_TF = cand_TF,
               cand_reciprocal = cand_reciprocal,
               relaxed_act = relaxed_act, relaxed_rep = relaxed_rep),
        file.path(rdsdir, "05_targets.rds"))

## 方法说明落盘 (论文 M&M 可直接引用)
methods_txt <- sprintf(
  paste0(
    "Candidate target identification (TF)\n",
    "Strategy: reciprocal-regulation screen using overexpression (OE) vs CK and RNAi knockdown vs CK.\n",
    "  Activated target: significant UP in OE_vs_CK AND significant DOWN in RNAi_vs_CK.\n",
    "  Repressed target: significant DOWN in OE_vs_CK AND significant UP in RNAi_vs_CK.\n",
    "Significance: padj < %.3f and |log2FoldChange| >= %.1f (shrunken LFC from DESeq2).\n",
    "This conservative double-threshold controls false positives but may exclude genes with a strong trend on only one side (see exploratory list and candidate list below).\n",
    "Custom score (heuristic, NOT a published statistic): score = log2FC(OE/CK) - log2FC(RNAi/CK).\n",
    "A true target is pushed in opposite directions by the two perturbations, so the difference equals\n",
    "the sum of the two effects. CAVEAT: the score is dominated by whichever side is larger and therefore\n",
    "ranks two biologically different classes together - genes changed only in OE (transgene/insertion\n",
    "or OE-specific induction) and genes changed only in RNAi (knockdown-specific candidate targets).\n",
    "The score alone therefore carries no evidence: tables/candidate_reciprocal_score.csv reports every\n",
    "gene with its score PLUS baseMean, log2FC/padj on both sides and an 'evidence' class\n",
    "(OE_only / RNAi_only / both_signif_opposite / both_signif_same_direction / neither_signif).\n",
    "tables/candidate_TF_locus.csv contains only the OE_only subset (significant change in OE and no\n",
    "significant change in RNAi), i.e. the class consistent with an over-expression locus/transgene;\n",
    "candidate identity still requires BLAST confirmation.\n"),
  padj_thresh, lfc_thresh)
writeLines(methods_txt, file.path(outdir, "tables", "targets_methods.txt"))

cat(sprintf("[%s] 严格靶标: 激活=%d, 抑制=%d; 探索性(宽松): 激活=%d, 抑制=%d (padj<%.3f, |LFC|>=%.1f, top_n=%d)\n",
            Sys.time(), length(act), length(rep),
            length(relaxed_act), length(relaxed_rep),
            padj_thresh, lfc_thresh, top_n))
cat(sprintf("  候选表: reciprocal_score=%d 行 (含 evidence 列), TF_locus(OE_only)=%d 行\n",
            nrow(cand_reciprocal), nrow(cand_TF)))
