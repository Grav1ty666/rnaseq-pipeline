# ============================================================
# 06_go.R  —  GO 注释 (必要步骤, 置于 07 富集之前)
#   产出 "基因->GO" 映射 (06_go.rds) 供 07_enrichment 消费。
#   两种来源:
#     * 默认: UniProt REST 在线抓取(按 gene_description 的 accession), 缓存 06_go_annotations.rds
#     * 离线: 设 GO_ANNOT_FILE=<文件> 直接读本地注释(如目标物种原生, 见 build_gene_go.R)
#   说明: GO 注释为同源转移(UniProt)或本地提供; 术语名/命名空间尽量保留。
#   失败处理: 解析/抓取失败一律 stop() (不再静默产出空表), 避免"注释缺失"被下游
#             误读成"无显著富集"; 06_go.rds 内含 status/source 供 07/09/10 核对。
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

f <- readRDS(file.path(rdsdir, "02_filtered.rds")); anno <- f$anno

## ---------- 本地 GO 注释解析(离线): GO_ANNOT_FILE 指定 TSV/CSV ----------
##   自动识别 gene/GO 列, 支持一格多个 GO(逗号/分号/竖线, 或 "term [GO:xxxxxxx]")。
parse_go_annot_file <- function(path) {
  ## 面向"干净的 gene->GO 文件"(带表头, 如 08/build_gene_go 的输出或用户自备);
  ## read.delim 默认 header=TRUE, 按表头名取 gene/GO/term/namespace 列。
  ## (eggNOG 原始注释不走这里, 而由 build_gene_go.R 先转成干净文件。)
  hdr <- readLines(path, n = 3, warn = FALSE)
  sep <- if (any(grepl("\t", hdr))) "\t" else ","
  df <- if (sep == "\t")
          read.delim(path, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE,
                     quote = "", na.strings = c("NA", ""))
        else
          read.csv(path, stringsAsFactors = FALSE, check.names = FALSE,
                   quote = "", na.strings = c("NA", ""))
  nm <- tolower(gsub("[^a-z0-9]+", "_", colnames(df)))
  pick <- function(pats) { for (p in pats) { i <- grep(p, nm); if (length(i)) return(i[1]) }; NA_integer_ }
  gc <- pick(c("^gene_id$", "^geneid$", "^gene$", "^gene_name$", "^genename$",
               "^protein_id$", "^proteinid$", "^transcript_id$", "^query$",
               "^query_id$", "^locus_tag$", "^locus$", "^id$"))
  oc <- pick(c("^go_id$", "^go_ids$", "^go$", "^gos$", "^go_term$", "^go_terms$",
               "^go_biological_process$", "^go_molecular_function$",
               "^go_cellular_component$", "^go_bp$", "^go_mf$", "^go_cc$"))
  if (is.na(gc) || is.na(oc)) {
    ## 退而求其次: 第1列为基因, GO 列取内容含 GO:xxxxxxx 的第一列
    gc <- 1
    oc <- NA_integer_
    for (j in seq_len(ncol(df))) { v <- head(as.character(df[[j]]), 500)
      if (any(grepl("GO:\\d{7}", v[!is.na(v)]))) { oc <- j; break } }
    if (is.na(oc))
      stop("本地 GO 文件无法识别 gene / GO 列。当前列: ", paste(colnames(df), collapse = ", "))
  }
  tc <- pick(c("^go_term$", "^go_name$", "^term$", "^description$", "^go_description$"))
  nc <- pick(c("^namespace$", "^ontology$", "^aspect$", "^type$", "^category$", "^go_domain$"))
  ## 向量化展开: 逐行构造 data.frame 在百万级注释上要跑数分钟, 这里一次性展开成列
  gene <- trimws(as.character(df[[gc]]))
  cell <- as.character(df[[oc]])
  ids_list <- regmatches(cell, gregexpr("GO:\\d{7}", cell))
  keep <- !is.na(cell) & nzchar(trimws(cell)) & lengths(ids_list) > 0
  ids_list <- lapply(ids_list[keep], unique)
  n_go <- lengths(ids_list)
  term_cell <- if (!is.na(tc)) trimws(as.character(df[[tc]])[keep]) else
               rep(NA_character_, sum(keep))
  nsp_cell  <- if (!is.na(nc)) as.character(df[[nc]])[keep] else
               rep(NA_character_, sum(keep))
  ## 一行多个 GO 时不使用该行的 term 名称 (原逻辑), 由下方以 GO id 兜底
  g <- unique(data.frame(
    gene_id   = rep(gene[keep], n_go),
    go_id     = unlist(ids_list, use.names = FALSE),
    go_term   = ifelse(n_go == 1, term_cell, NA_character_),
    namespace = rep(nsp_cell, n_go),
    stringsAsFactors = FALSE))
  miss <- is.na(g$go_term) | !nzchar(trimws(g$go_term)); g$go_term[miss] <- g$go_id[miss]
  ns <- gsub("[^a-z]", "", tolower(g$namespace))
  g$namespace <- ifelse(grepl("biol|process|^bp$", ns), "BP",
                 ifelse(grepl("cell|component|^cc$", ns), "CC",
                 ifelse(grepl("molecular|function|^mf$", ns), "MF", "NA")))
  n_match <- length(intersect(g$gene_id, anno$gene_id))
  cat(sprintf("[GO] 本地注释(离线): %d 条 gene-term, 唯一 GO %d; 与本研究基因匹配 %d/%d\n",
              nrow(g), length(unique(g$go_id)), n_match, length(anno$gene_id)))
  if (nrow(g) == 0) stop("[GO] 本地注释文件未解析出任何 gene-term, 请检查格式/GO ID 列。")
  if (n_match < 0.1 * length(anno$gene_id))
    warning("[GO] 本地注释与本研究 gene_id 匹配率低(<10%); 请确认基因 ID 与 count 矩阵 gene_id 一致。")
  list(go_term2gene = g, go_name = setNames(g$go_term, g$go_id), gene_acc = NULL,
       n_uniq_acc = length(unique(g$gene_id)), kegg_term2gene = NULL, kegg_available = FALSE,
       source = "local", source_file = path, status = "ok")
}

## ---------- 获取 GO 映射: 优先本地文件, 否则 UniProt 在线抓取(缓存) ----------
get_go_annotations <- function() {
  ## 默认: 优先本地 GO 注释文件(若存在), 否则回退 UniProt 抓取;
  ##       设 GO_ANNOT_FILE 可显式指定其它本地 GO 文件。
  go_file <- Sys.getenv("GO_ANNOT_FILE", "")
  if (!nzchar(go_file)) {
    local_go <- Sys.getenv("LOCAL_GO_FILE", "~/eggnog_data/go_by_gene.tsv")
    if (file.exists(local_go)) {
      cat(sprintf("[GO] 默认使用本地 GO 文件: %s\n", local_go)); go_file <- local_go
    } else {
      cat(sprintf("[GO] 未设 GO_ANNOT_FILE 且本地 GO 文件不存在(%s), 回退 UniProt 抓取\n", local_go))
    }
  }
  if (nzchar(go_file)) {
    if (!file.exists(go_file)) stop("GO 文件不存在: ", go_file)
    return(parse_go_annot_file(go_file))
  }
  out_rds <- file.path(rdsdir, "06_go_annotations.rds")
  if (file.exists(out_rds) && Sys.getenv("FORCE_GO") != "1") {
    cat(sprintf("[GO] 复用已有 GO 映射 (%s)\n", out_rds)); return(readRDS(out_rds))
  }
  parse_acc <- function(desc) {
    if (is.na(desc) || desc == "") return(character(0))
    m <- regmatches(desc, gregexpr("(?:sp|tr)\\|([A-Z0-9]+)\\|", desc, perl = TRUE))[[1]]
    if (length(m) == 0) return(character(0))
    unique(sapply(strsplit(m, "\\|"), `[`, 2))
  }
  gene_acc <- setNames(lapply(anno$gene_description, parse_acc), anno$gene_id)
  gene_acc <- gene_acc[sapply(gene_acc, length) > 0]
  uniq_acc <- unique(unlist(gene_acc))
  cat(sprintf("[GO] 保留基因 %d, 带 UniProt 命中 %d, 唯一 accession %d\n",
              nrow(anno), length(gene_acc), length(uniq_acc)))
  curl_tsv <- function(url, timeout = 45) {
    tf <- tempfile(fileext = ".tsv")
    st <- system2("curl", c("-s", "--max-time", as.character(timeout), "-o", tf, shQuote(url)),
                 stdout = FALSE, stderr = FALSE)
    if (st != 0 || !file.exists(tf) || file.size(tf) == 0) { unlink(tf); return(NULL) }
    d <- tryCatch(read.delim(tf, stringsAsFactors = FALSE, check.names = FALSE), error = function(e) NULL)
    unlink(tf); d
  }
  parse_go_cell <- function(cell) {
    if (is.na(cell) || cell == "") return(data.frame(go_id = character(0), go_term = character(0)))
    terms <- regmatches(cell, gregexpr("(.+?)\\s+\\[(GO:\\d+)\\]", cell))[[1]]
    if (length(terms) == 0) return(data.frame(go_id = character(0), go_term = character(0)))
    data.frame(go_id  = sub("^.*\\[(GO:\\d+)\\]$", "\\1", terms),
               go_term = trimws(gsub("^[;[:space:]]+", "", sub("\\s+\\[GO:\\d+\\]$", "", terms))),
               stringsAsFactors = FALSE)
  }
  go_cols <- c("Gene Ontology (biological process)" = "BP",
               "Gene Ontology (cellular component)" = "CC",
               "Gene Ontology (molecular function)" = "MF")
  acc_to_genes <- list()
  for (g in names(gene_acc)) for (a in gene_acc[[g]]) acc_to_genes[[a]] <- c(acc_to_genes[[a]], g)
  chunk <- 100; rows <- list(); failed <- 0
  for (i in seq(1, length(uniq_acc), chunk)) {
    accs <- uniq_acc[i:min(i + chunk - 1, length(uniq_acc))]
    qry  <- paste0("accession:(", paste(accs, collapse = " OR "), ")")
    url  <- paste0("https://rest.uniprot.org/uniprotkb/search?query=", gsub(" ", "%20", qry),
                   "&fields=accession,go_p,go_c,go_f&format=tsv&size=500")
    d <- curl_tsv(url)
    if (is.null(d) || nrow(d) == 0 || !"Entry" %in% colnames(d)) { failed <- failed + 1; Sys.sleep(0.2); next }
    d$Entry <- sub("\\..*$", "", as.character(d$Entry))
    for (col in names(go_cols)) {
      ns <- go_cols[[col]]
      for (j in seq_len(nrow(d))) {
        a <- d$Entry[j]; gs <- acc_to_genes[[a]]
        if (length(gs) == 0) next
        p <- parse_go_cell(d[j, col])
        if (nrow(p)) rows <- c(rows, list(data.frame(
          gene_id = rep(gs, each = nrow(p)), go_id = rep(p$go_id, length(gs)),
          go_term = rep(p$go_term, length(gs)), namespace = ns, stringsAsFactors = FALSE)))
      }
    }
    Sys.sleep(0.15)
  }
  go_term2gene <- if (length(rows)) unique(do.call(rbind, rows)) else
    data.frame(gene_id = character(0), go_id = character(0), go_term = character(0),
               namespace = character(0), stringsAsFactors = FALSE)
  ## 抓取失败必须显式报错: 空 GO 表会被 07 跳过、进而被报告误读成"无显著富集"
  if (nrow(go_term2gene) == 0)
    stop("[GO] UniProt 抓取未获得任何 gene-term (网络不可达或 accession 无命中)。\n",
         "  请提供本地注释后重跑 06: GO_ANNOT_FILE=<gene->GO 文件> Rscript ", "Pipeline/06_go.R",
         "\n  (不要据此认为'无富集'。)")
  go_name <- setNames(go_term2gene$go_term, go_term2gene$go_id)
  cat(sprintf("[GO] 映射: %d 条 (gene-term), 唯一 GO %d, 抓取失败块 %d\n",
              nrow(go_term2gene), length(unique(go_term2gene$go_id)), failed))
  kegg_ok <- function() {
    tf <- tempfile()
    st <- system2("curl", c("-s", "--max-time", "15", "-o", tf, "https://rest.kegg.jp/conv/ko/P04637"),
                 stdout = FALSE, stderr = FALSE)
    ok <- st == 0 && file.exists(tf) && file.size(tf) > 0; unlink(tf); ok
  }
  kegg_term2gene <- NULL
  if (kegg_ok()) { cat("[GO] 检测到 KEGG 可达\n"); warning("KEGG 抓取接口已预留, 当前未实现完整映射") }
  else cat("[GO] 警告: 本环境无法访问 kegg.jp, 跳过 KEGG 富集\n")
  res <- list(go_term2gene = go_term2gene, go_name = go_name, gene_acc = gene_acc,
              n_uniq_acc = length(uniq_acc), kegg_term2gene = kegg_term2gene,
              kegg_available = !is.null(kegg_term2gene),
              source = "uniprot", source_file = NA_character_, status = "ok")
  saveRDS(res, file = out_rds)
  cat(sprintf("[GO] GO 注释映射已保存 (%s)\n", out_rds))
  res
}

## 失败即中断: 否则空的 06_go.rds 会让 07 跳过 GO 富集, 而报告里"0 项显著"会被
## 误读为"看过了、没有富集"(这正是之前踩过的坑)。
go_ann <- get_go_annotations()
if (is.null(go_ann$status)) go_ann$status <- "ok"
saveRDS(go_ann, file.path(rdsdir, "06_go.rds"))

## 注释来源留痕 (路径/大小/md5/时间/规模): 论文方法学需要能追溯到具体注释文件版本
prov <- c(
  sprintf("GO annotation source: %s (%s)", go_ann$source, Sys.time()),
  if (!is.null(go_ann$source_file) && !is.na(go_ann$source_file) && file.exists(go_ann$source_file))
    c(sprintf("  file      : %s", go_ann$source_file),
      sprintf("  size      : %s bytes", format(file.size(go_ann$source_file), big.mark = ",")),
      sprintf("  md5       : %s", unname(tools::md5sum(go_ann$source_file))),
      sprintf("  mtime     : %s", format(file.mtime(go_ann$source_file))))
  else "  file      : (UniProt REST 抓取, 见 06_go_annotations.rds; 未记录抓取日期/返回批次)",
  sprintf("  gene-term : %d 条; 基因 %d 个; 唯一 GO %d 个",
          nrow(go_ann$go_term2gene), length(unique(go_ann$go_term2gene$gene_id)),
          length(unique(go_ann$go_term2gene$go_id))),
  sprintf("  R/DESeq2  : R %s.%s; DESeq2 %s", R.version$major, R.version$minor,
          as.character(utils::packageVersion("DESeq2"))),
  "  提醒: 该注释文件位于工程目录之外, 仅凭本仓库无法复现 GO 结果; 建议随补充材料一并归档")
writeLines(prov, file.path(outdir, "tables", "GO_source_provenance.txt"))

cat(sprintf("[%s] GO 注释步骤完成 (source=%s): %d 条 gene-term, %d 唯一 GO, 覆盖 %d 个基因 (已存 06_go.rds, 供 07 富集使用)\n",
            Sys.time(), go_ann$source, nrow(go_ann$go_term2gene),
            length(unique(go_ann$go_term2gene$go_id)),
            length(unique(go_ann$go_term2gene$gene_id))))
