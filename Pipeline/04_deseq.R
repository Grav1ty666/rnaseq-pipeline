# ============================================================
# 04_deseq.R  —  DESeq2 差异分析 (三组两两对比) + DEG 列表
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

## ---------- 参数 (可经环境变量覆盖, 调参无需改源码) ----------
padj_thresh <- as.numeric(Sys.getenv("DEG_PADJ", "0.05"))
lfc_thresh  <- as.numeric(Sys.getenv("DEG_LFC",  "1"))
design_str  <- Sys.getenv("DESEQ_DESIGN", "~ condition")  # 可选: "~ batch + condition"

f <- readRDS(file.path(rdsdir, "02_filtered.rds"))
count_mat <- f$count_mat
anno      <- f$anno

## ---------- 构建 colData ----------
coldat <- data.frame(row.names = sample_cols, condition = cond)

## 5) coldat 行名 与 count_mat 列名一致性校验
if (!all(colnames(count_mat) %in% rownames(coldat)))
  stop("count_mat 列名与 sample_cols 不一致, 请检查输入/上游 01。")
if (!identical(colnames(count_mat), rownames(coldat)))
  warning("count_mat 列序与 sample_cols 不同; DESeq2 按名称匹配, 请确认分组正确。")

## 1) 可选批次校正: 经 DESEQ_BATCH 传入与 sample_cols 同序的批次向量
batch_env <- Sys.getenv("DESEQ_BATCH", "")
if (batch_env != "") {
  batch_vec <- strsplit(batch_env, ",")[[1]]
  if (length(batch_vec) != length(sample_cols))
    stop("DESEQ_BATCH 长度须等于样本数 (", length(sample_cols), ")")
  coldat$batch <- factor(batch_vec)
  design_str <- "~ batch + condition"
  cat("注意: 已启用批次校正, 设计 = ~ batch + condition\n")
}

## 2) CK 仅 2 重复: 离散度估计不稳, 提醒
if (sum(cond == "CK") < 3)
  warning("CK 仅 ", sum(cond == "CK"),
          " 个生物学重复: DESeq2 离散度估计不稳, padj 波动大, 假阳性风险升高, 解读需谨慎。")
## 1) 未校正批次的提醒
if (!grepl("batch", design_str))
  cat("提示: 当前设计 = ", design_str,
      "; 若测序/建库存在批次差异, 结果会混杂批次影响。",
      "可用 DESEQ_BATCH 提供批次变量启用 ~ batch + condition 校正。\n", sep = "")

design_f <- as.formula(design_str)
dds <- DESeqDataSetFromMatrix(countData = count_mat,
                              colData   = coldat,
                              design    = design_f)
dds <- DESeq(dds)
cat(sprintf("[%s] DESeq2 完成\n", Sys.time()))

## 6) lfcShrink: LFC/lfcSE 用收缩后验值, pvalue/stat/padj 仍来自未收缩的 Wald 检验
##    (与 DESeq2 官方 lfcShrink 输出口径一致: 收缩不改变 p 值)
##    apeglm 只支持 coef=, 不支持 contrast=, 因此:
##      * OE_vs_CK / RNAi_vs_CK : coef 收缩 (apeglm), LFC 与 lfcSE 均为后验值
##      * OE_vs_RNAi            : 默认用两 coef 后验均值之差 (posterior mean 的线性性)。
##        该"差值"没有可用的后验 SE, 故只替换 LFC, lfcSE 保留 MLE 值 ->
##        **lfcSE/stat/pvalue 三者与收缩后的 LFC 不可互相换算** (已写入 DEG_methods.txt)。
##        若需要该对比也具备自洽的 LFC/lfcSE, 设 DESEQ_SHRINK=normal
##        (改用 contrast 式 normal 收缩; 缺点: normal 先验偏差大于 apeglm)。
shrink_env  <- Sys.getenv("DESEQ_SHRINK", "apeglm")   # apeglm | normal
shrink_type <- if (requireNamespace("apeglm", quietly = TRUE) && shrink_env != "normal")
                 "apeglm" else "normal"
shrink_methods <- character(0)   # 逐对比记录实际所用收缩方式

shrink_for <- function(contrast) {
  lvl2 <- contrast[2]; lvl1 <- contrast[3]
  nm   <- paste(lvl2, "vs", lvl1)
  if (shrink_type == "apeglm") {
    a <- lfcShrink(dds, coef = paste0("condition_", lvl2, "_vs_CK"), type = "apeglm")
    if (lvl1 == "CK") {
      shrink_methods[[nm]] <<- "apeglm (coef 收缩; log2FoldChange 与 lfcSE 均为后验值)"
      return(list(log2FoldChange = a$log2FoldChange, lfcSE = a$lfcSE))
    }
    b <- lfcShrink(dds, coef = paste0("condition_", lvl1, "_vs_CK"), type = "apeglm")
    shrink_methods[[nm]] <<-
      "apeglm (两 coef 后验均值之差; 该对比无后验 SE, lfcSE/stat/pvalue 仍为 MLE 值)"
    return(list(log2FoldChange = a$log2FoldChange - b$log2FoldChange, lfcSE = NULL))
  }
  out <- lfcShrink(dds, contrast = contrast, type = "normal")
  shrink_methods[[nm]] <<- "normal (contrast 收缩; log2FoldChange 与 lfcSE 均为后验值)"
  list(log2FoldChange = out$log2FoldChange, lfcSE = out$lfcSE)
}
safe_shrink <- function(contrast) {
  nm <- paste(contrast[2], "vs", contrast[3])
  tryCatch(shrink_for(contrast),
           error = function(e) {
             warning("lfcShrink(", shrink_type, ") 对 ", nm, " 失败, 回退 normal: ",
                     conditionMessage(e))
             shrink_type <<- "normal"
             out <- lfcShrink(dds, contrast = contrast, type = "normal")
             shrink_methods[[nm]] <<- "normal (回退: 初次收缩失败)"
             list(log2FoldChange = out$log2FoldChange, lfcSE = out$lfcSE)
           })
}
apply_shrink <- function(res, contrast) {
  s <- safe_shrink(contrast)
  res$log2FoldChange <- s$log2FoldChange
  if (!is.null(s$lfcSE)) res$lfcSE <- s$lfcSE   # 仅当收缩给出后验 SE 时一并替换
  res
}
res_OE_CK   <- apply_shrink(results(dds, contrast = c("condition", "OE",   "CK"), alpha = padj_thresh),
                            c("condition", "OE",   "CK"))
res_RNAi_CK <- apply_shrink(results(dds, contrast = c("condition", "RNAi", "CK"), alpha = padj_thresh),
                            c("condition", "RNAi", "CK"))
res_OE_RNAi <- apply_shrink(results(dds, contrast = c("condition", "OE",   "RNAi"), alpha = padj_thresh),
                            c("condition", "OE",   "RNAi"))

## ---------- 整理为带注释的表格 ----------
## 4) 注释对齐前防御: anno$gene_id 必须唯一 (否则 merge 会膨胀行数)
if (any(duplicated(anno$gene_id)))
  stop("anno$gene_id 存在重复, merge 注释会错乱; 请先修复上游 01/02。")

build_tab <- function(res, name) {
  t <- as.data.frame(res)
  t$gene_id <- rownames(t)
  ## 防 merge 列名冲突: anno 若含 DESeq 结果同名列会生成 .x/.y 后缀, 令下游静默失效
  reserved <- c("log2FoldChange","padj","baseMean","pvalue","lfcSE","stat")
  clash <- intersect(reserved, setdiff(colnames(anno), "gene_id"))
  if (length(clash) > 0)
    stop("注释表含与 DESeq 结果冲突的列名: ", paste(clash, collapse = ", "),
         "; 请在输入前重命名, 否则 merge 会产生 .x/.y 后缀导致下游取列失败。")
  t <- merge(t, anno, by = "gene_id", all.x = TRUE)   # merge 会按 gene_id 重排, 下面再排序
  t <- t[order(t$padj, -abs(t$log2FoldChange), na.last = TRUE), ]
  write.csv(t, file.path(outdir, "tables", sprintf("DEG_%s.csv", name)),
            row.names = FALSE)
  t
}
t_OE_CK   <- build_tab(res_OE_CK,   "OE_vs_CK")
t_RNAi_CK <- build_tab(res_RNAi_CK, "RNAi_vs_CK")
t_OE_RNAi <- build_tab(res_OE_RNAi, "OE_vs_RNAi")

## 3) DEG 判定 (参数化阈值; 阈值常规标准, 方法部分已声明)
deg <- function(t) t$gene_id[!is.na(t$padj) & t$padj < padj_thresh &
                             !is.na(t$log2FoldChange) & abs(t$log2FoldChange) >= lfc_thresh]
deg_OE   <- deg(t_OE_CK);   n_OE   <- length(deg_OE)
deg_RNAi <- deg(t_RNAi_CK); n_RNAi <- length(deg_RNAi)
deg_OR   <- deg(t_OE_RNAi); n_OR   <- length(deg_OR)

## ---------- 保存 (含原始 res 对象, 便于改阈值/重提取不必重跑 DESeq) ----------
saveRDS(list(dds = dds,
             t_OE_CK = t_OE_CK, t_RNAi_CK = t_RNAi_CK, t_OE_RNAi = t_OE_RNAi,
             res_OE_CK = res_OE_CK, res_RNAi_CK = res_RNAi_CK, res_OE_RNAi = res_OE_RNAi,
             deg_OE = deg_OE, deg_RNAi = deg_RNAi, deg_OR = deg_OR,
             n_OE = n_OE, n_RNAi = n_RNAi, n_OR = n_OR,
             params = list(padj_thresh = padj_thresh, lfc_thresh = lfc_thresh,
                           design = design_str, shrink_type = shrink_type,
                           shrink_methods = shrink_methods)),
        file.path(rdsdir, "04_deseq.rds"))

## ---------- 方法说明落盘 (论文 M&M 可直接引用) ----------
methods_txt <- c(
  "Differential expression analysis",
  sprintf("Method: DESeq2 (Love et al., 2014), design = %s.", design_str),
  sprintf("Significance: padj < %.3f and |log2FoldChange| >= %.1f (DESeq2 independent filtering on).",
          padj_thresh, lfc_thresh),
  sprintf("log2 fold changes were shrunk with lfcShrink (overall type = %s; apeglm: Zhu et al. 2019).",
          shrink_type),
  "Shrinkage actually applied per contrast:",
  sprintf("  %s: %s", names(shrink_methods), unlist(shrink_methods, use.names = FALSE)),
  "IMPORTANT: log2FoldChange and lfcSE in DEG_*.csv are posterior (shrunken) values, while",
  "stat/pvalue/padj come from the unshrunken Wald test (this is the standard DESeq2 lfcShrink",
  "output convention). For contrasts labelled \"no posterior SE\", lfcSE/stat/pvalue are MLE-based",
  "and must not be combined with the shrunken LFC (e.g. do not recompute z = log2FoldChange/lfcSE).",
  "Alternative: set DESEQ_SHRINK=normal for contrast-based normal shrinkage, which returns a",
  "consistent log2FoldChange + lfcSE for every contrast (at the cost of a slightly more biased LFC).",
  "type='ashr' would also give posterior SEs for custom contrasts but requires the ashr package,",
  "which is not installed in this environment.")
writeLines(methods_txt, file.path(outdir, "tables", "DEG_methods.txt"))

cat(sprintf("[%s] DEG: OE_vs_CK=%d, RNAi_vs_CK=%d, OE_vs_RNAi=%d (padj<%.3f, |LFC|>=%.1f, shrink=%s)\n",
            Sys.time(), n_OE, n_RNAi, n_OR, padj_thresh, lfc_thresh, shrink_type))
cat("  收缩方式: ", paste(sprintf("%s=%s", names(shrink_methods), unlist(shrink_methods, use.names = FALSE)),
                          collapse = "; "), "\n", sep = "")
