# ============================================================
# 01_read_data.R  —  读取原始 count 矩阵 + 注释
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

## 1) 文件存在性前置判断 (错误信息更直观)
if (!file.exists(countf))
  stop("count 文件不存在: ", countf,
       "\n请检查 RNASEQ_COUNT 环境变量或文件路径是否正确。")

df <- read.delim(countf, stringsAsFactors = FALSE, check.names = FALSE)

## 2) gene_id 列校验: 改名会静默出错; 重复 ID 会破坏 DESeq2 唯一行名
geneid_col <- "gene_id"
if (!(geneid_col %in% colnames(df)))
  stop("count 文件缺少 '", geneid_col, "' 列; 现有前几列: ",
       paste(colnames(df)[1:min(10, length(colnames(df)))], collapse = ", "), " ...")
if (any(duplicated(df[[geneid_col]])))
  stop("count 文件存在重复的 ", geneid_col, " (共 ",
       sum(duplicated(df[[geneid_col]])), " 个); DESeq2 要求唯一行名, 请先合并/去重。")

## 3) 样本列存在性校验 (df[, sample_cols] 会按 sample_cols 顺序重排, 故只校验存在性)
miss <- setdiff(sample_cols, colnames(df))
if (length(miss) > 0)
  stop("count 文件缺少样本列: ", paste(miss, collapse = ", "),
       "\n  现有列: ", paste(colnames(df), collapse = ", "))

count_mat <- as.matrix(df[, sample_cols])
rownames(count_mat) <- df[[geneid_col]]

## 4) 非整数计数告警 (storage.mode 强制转 integer 会静默截断)
if (any(count_mat != round(count_mat)))
  warning("count 矩阵存在非整数值, 将被强制截断为整数 (可能影响下游统计准确性)。")

storage.mode(count_mat) <- "integer"

## 保留完整注释信息(后续过滤步骤再做对齐)
anno_full <- df

saveRDS(list(count_mat = count_mat, anno_full = anno_full),
        file.path(rdsdir, "01_raw.rds"))

cat(sprintf("[%s] 原始基因数: %d\n", Sys.time(), nrow(count_mat)))
