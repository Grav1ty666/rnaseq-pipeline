# ============================================================
# 02_filter.R  —  低表达过滤 (>=min_count 的样本数 >= min_samples) + 注释对齐
# 过滤规则: 至少 min_samples 个样本中 read 计数 >= min_count (转录组常规策略)
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

## 1) 过滤参数 (可经环境变量覆盖, 调参无需改源码)
min_count   <- as.integer(Sys.getenv("FILTER_MIN_COUNT",   "10"))
min_samples <- as.integer(Sys.getenv("FILTER_MIN_SAMPLES", "2"))
family_col  <- Sys.getenv("FAMILY_COL", "Family")

raw <- readRDS(file.path(rdsdir, "01_raw.rds"))
count_mat <- raw$count_mat
anno_full <- raw$anno_full

## 2) 防御: anno_full$gene_id 必须唯一 (01 已校验, 此处再确认, 否则 match 对齐错乱)
if (any(duplicated(anno_full$gene_id)))
  stop("anno_full$gene_id 存在重复, 注释对齐将错乱; 请先修复输入/上游 01。")

## 3) 低表达过滤
keep <- rowSums(count_mat >= min_count) >= min_samples
count_mat <- count_mat[keep, , drop = FALSE]

## 4) 注释按保留基因对齐 (match 1:1, 并检查是否有 NA)
m <- match(rownames(count_mat), anno_full$gene_id)
if (any(is.na(m)))
  stop("注释对齐失败: ", sum(is.na(m)),
       " 个基因在 anno_full 中无对应 gene_id (count 与注释行数/ID 不一致)。")
anno <- anno_full[m, , drop = FALSE]
anno$gene_id <- rownames(count_mat)

## 5) Family 列存在性校验 + 规范化到统一列名 'Family'
##    下游 05/07/08/09 一律按 'Family' 读取, 故 FAMILY_COL 只能在这里生效一次:
##    被选中的列会复制为 'Family', 避免"只有 02 认识 FAMILY_COL、下游仍读 Family"
##    导致的隐性失效(下游会读到未清洗的 NA 甚至报错)。
if (!(family_col %in% colnames(anno)))
  stop("注释表缺少 '", family_col, "' 列; 现有列: ",
       paste(colnames(anno), collapse = ", "))
if (!identical(family_col, "Family")) {
  cat(sprintf("[02] FAMILY_COL=%s -> 规范化列名 'Family' (供下游 05/07/08/09 统一读取)\n",
              family_col))
  anno$Family <- anno[[family_col]]
}
anno$Family[is.na(anno$Family)] <- "-"

## 5b) 下游依赖列校验 (gene_name/gene_description 被 04/05/07 使用, 缺失即明确报错)
for (need in c("gene_name", "gene_description"))
  if (!(need %in% colnames(anno)))
    stop("注释表缺少 '", need, "' 列 (下游 04/05/07 依赖); 现有列: ",
         paste(colnames(anno), collapse = ", "))

## 6) 对齐后一致性校验
if (!identical(anno$gene_id, rownames(count_mat)))
  stop("对齐后 anno$gene_id 与 count_mat 行名不一致, 请检查。")

## params 记录本次过滤/注释口径, 供下游与归档核对 (下游仍只按名称取 count_mat/anno)
saveRDS(list(count_mat = count_mat, anno = anno,
             params = list(family_col = family_col)),
        file.path(rdsdir, "02_filtered.rds"))

## 7) 过滤统计落盘 (参数 + 计数, 便于复现归档)
stats <- data.frame(
  parameter = c("min_count", "min_samples", "n_before", "n_kept", "n_removed"),
  value     = c(min_count, min_samples,
               nrow(raw$count_mat), nrow(count_mat), sum(!keep)),
  stringsAsFactors = FALSE
)
write.csv(stats, file.path(outdir, "tables", "filter_stats.csv"), row.names = FALSE)

## 落盘预处理方法说明 (含文献), 覆盖 01(导入) 与 02(过滤)
preproc <- c(
  "Preprocessing methods (TF RNA-seq): read import (01) and low-count filter (02)",
  "==================================================================================",
  "1) Read import (01_read_data.R)",
  "   - Input: per-gene raw read-count matrix (gene_count.xls, tab-delimited): 'gene_id' column",
  "     + one integer count column per sample.",
  "   - Raw counts assumed from a standard RNA-seq quantification pipeline: alignment with STAR",
  "     [Dobin et al. 2013] or HISAT2 [Kim et al. 2015], then read counting with featureCounts",
  "     [Liao, Smyth & Shi 2014] or HTSeq [Anders et al. 2015]. This pipeline consumes those",
  "     counts; it does not re-quantify.",
  "   - Checks: gene_id present & unique (DESeq2 requires unique row names); sample columns present;",
  "     non-integer counts flagged.",
  "",
  "2) Low-count filtering (02_filter.R)",
  "   - Rule: keep a gene if count >= min_count (default 10) in >= min_samples (default 2) samples.",
  "     Standard low-count removal before DE to reduce noise and aid dispersion estimation",
  "     [Chen et al. 2016; Law et al. 2014].",
  "   - Count matrix re-aligned 1:1 to annotation by gene_id; TF 'Family' and downstream columns",
  "     (gene_name / gene_description) validated. Stats in tables/filter_stats.csv.",
  "   - The TF family column can be pointed to another annotation column with FAMILY_COL; it is then",
  "     normalised to the canonical name 'Family' so that every downstream step reads the same column.",
  "",
  "References:",
  "  Dobin A, et al. STAR: ultrafast universal RNA-seq aligner. Bioinformatics. 2013;29(1):15-21.",
  "    DOI:10.1093/bioinformatics/bts635 | PMID:23104886",
  "  Kim D, Langmead B, Salzberg SL. HISAT: a fast spliced aligner with low memory requirements.",
  "    Nat Methods. 2015;12(4):357-360. DOI:10.1038/nmeth.3317 | PMID:25751142",
  "  Liao Y, Smyth GK, Shi W. featureCounts. Bioinformatics. 2014;30(7):923-930.",
  "    DOI:10.1093/bioinformatics/btt656 | PMID:24227677",
  "  Anders S, Pyl PT, Huber W. HTSeq. Bioinformatics. 2015;31(2):166-169. DOI:10.1093/bioinformatics/btu638 | PMID:25260700",
  "  Chen Y, Lun ATL, Smyth GK. Differential expression analysis of complex RNA-seq experiments using",
  "    edgeR. Stat Anal Genomics Data Protoc. 2016:51-74.",
  "  Law CW, Chen Y, Shi W, Smyth GK. voom. Genome Biol. 2014;15(2):R29. DOI:10.1186/gb-2014-15-2-r29 | PMID:24485249")
writeLines(preproc, file.path(outdir, "tables", "preprocessing_methods.txt"))

cat(sprintf("[%s] 过滤后保留基因数: %d (剔除 %d 个低/零表达; 阈值 >=%d reads in >=%d samples)\n",
            Sys.time(), nrow(count_mat), sum(!keep), min_count, min_samples))
