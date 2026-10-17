# ============================================================
# build_gene_go.R  —  用参考基因组注释生成本研究
#                   gene_id 口径的 GO 文件(离线, 替代 UniProt 同源转移)
# 前提(已验证): 本研究 gene-XXXXXXXX 与 gff gene ID 100% 一一对应,
#               蛋白 PROT##### -> 基因 可由 gff CDS(protein_id + locus_tag)提取。
# 输入(环境变量):
#   GFF_FILE     = ncbi_dataset/data/reference/genomic.gff
#   PROT_GO_FILE  = 蛋白(PROT)->GO 文件(eggNOG-mapper .emapper.annotations 或 InterProScan .tsv)
#   GO_OBO        = (可选) go-basic.obo, 用于补 GO term 名称与 BP/CC/MF 命名空间
# 输出:
#   results/tables/protein2gene.tsv  (PROT -> gene_id)
#   results/tables/go_by_gene.tsv    (gene_id, go_id, go_term, namespace)
# 用法:
#   # 1) 先对 protein.faa 跑 eggNOG-mapper(推荐) 或 InterProScan, 得到 PROT_GO_FILE
#   # 2) 生成本研究 GO 文件:
#   GFF_FILE=... PROT_GO_FILE=... Rscript Pipeline/build_gene_go.R
#   # 3) 离线富集:
#   GO_ANNOT_FILE=<go_by_gene.tsv> Rscript Pipeline/06_go.R     # GO 注释
#   Rscript Pipeline/07_enrichment.R                                  # 富集(依赖 06 产物)
# ============================================================
script_dir <- (function() {
  args <- commandArgs(trailingOnly = FALSE)
  fa <- grep("^--file=", args, value = TRUE)
  if (length(fa)) return(dirname(sub("^--file=", "", fa[1])))
  "."
})()
source(file.path(script_dir, "00_config.R"))

## 本研究基因(过滤后)来自 02; 用于统计映射覆盖率
if (!file.exists(file.path(rdsdir, "02_filtered.rds")))
  stop("缺少 02_filtered.rds, 请先运行 02_filter.R")
anno <- readRDS(file.path(rdsdir, "02_filtered.rds"))$anno

gff <- Sys.getenv("GFF_FILE",
        file.path(base, "data", "reference", "genomic.gff"))
if (!file.exists(gff)) stop("GFF_FILE 不存在: ", gff)

## ---- 1) 从 gff CDS 提取 protein_id(PROT) -> gene_id(gene-<locus_tag>) ----
xwalk <- file.path(outdir, "tables", "protein2gene.tsv")
cat(sprintf("[%s] 解析 gff CDS -> protein2gene (这步约十几秒)...\n", Sys.time()))
awk_cmd <- sprintf(
  "awk -F'\\t' '$3==\"CDS\"{pid=\"\";lt=\"\";if(match($9,/protein_id=[^;]+/))pid=substr($9,RSTART+11,RLENGTH-11);if(match($9,/locus_tag=[^;]+/))lt=substr($9,RSTART+10,RLENGTH-10);if(pid!=\"\"&&lt!=\"\")print pid\"\\tgene-\"lt}' %s | sort -u > %s",
  shQuote(gff), shQuote(xwalk))
st <- system(awk_cmd)
if (st != 0 || !file.exists(xwalk) || file.size(xwalk) == 0)
  stop("gff CDS 解析失败(protein2gene 为空)。请检查 gff 格式/列。")
xw <- read.delim(xwalk, header = FALSE, stringsAsFactors = FALSE,
                 col.names = c("protein_id", "gene_id"))
cat(sprintf("[annotation] protein2gene: %d 个蛋白 -> %d 个基因\n",
            nrow(xw), length(unique(xw$gene_id))))
n_in <- length(intersect(xw$gene_id, anno$gene_id))
cat(sprintf("[annotation] 其中与本研究 gene_id 匹配: %d 个基因\n", n_in))

## ---- 2) 若提供 PROT_GO_FILE, 转成本研究 gene_id 口径的 GO 文件 ----
prot_go <- Sys.getenv("PROT_GO_FILE",
                      "~/eggnog_data/emapper/out.emapper.annotations")
out_go  <- Sys.getenv("GO_OUT", "~/eggnog_data/go_by_gene.tsv")  # 持久位置, 作 06_go 默认 GO 源
if (!nzchar(prot_go)) {
  cat("\n[下一步] 尚未提供 PROT_GO_FILE。请对 protein.faa 跑:\n")
  cat("   eggNOG-mapper: emapper.py -i protein.faa -o reference --dmnd_db eggnog ... (得 out.emapper.annotations)\n")
  cat("   或 InterProScan: interproscan.sh -i protein.faa -f TSV (得 .tsv)\n")
  cat(sprintf("然后: GFF_FILE=%s PROT_GO_FILE=<你的文件> Rscript Pipeline/build_gene_go.R\n", gff))
  cat("(本次已生成 protein2gene 映射, 可先行核对)\n")
  quit(save = "no", status = 0)
}
if (!file.exists(prot_go)) stop("PROT_GO_FILE 不存在: ", prot_go)

## 读取蛋白->GO 文件(自动识别 蛋白ID列 / GO列 / 术语列; 兼容 eggNOG 与 InterProScan)
## 分隔符用"数据行(跳过 # 注释)"探测: eggNOG 首行是 ## 注释(无 tab), 不能拿它判断
raw_lines <- readLines(prot_go, n = 50, warn = FALSE)
data_lines <- raw_lines[!grepl("^#", raw_lines)]
sep <- if (any(grepl("\t", data_lines))) "\t" else ","
read_fun <- if (sep == "\t") read.delim else read.csv
pg <- read_fun(prot_go, header = FALSE, stringsAsFactors = FALSE, quote = "",
               na.strings = c("NA",""), fill = TRUE, comment.char = "#")  # comment.char="#" 跳过 ## 与 #query 头
rownames(pg) <- NULL
idc <- 1   # 第1列 = 蛋白/query ID (eggNOG 与 InterProScan 均如此)
## GO 列: 取内容含 GO:xxxxxxx 的第一列(仅扫前 500 行以加速)
goc <- NA_integer_
for (j in seq_len(ncol(pg))) { v <- head(as.character(pg[[j]]), 500); if (any(grepl("GO:\\d{7}", v[!is.na(v)]))) { goc <- j; break } }
if (is.na(goc)) stop("无法在 PROT_GO_FILE 中定位 GO 列(未找到 GO:xxxxxxx)。")

pwa2gene <- setNames(xw$gene_id, sub("\\..*$", "", xw$protein_id))  # 键去版本号
pwa2gene <- pwa2gene[!duplicated(names(pwa2gene))]
rows <- list()
pid_col <- trimws(as.character(pg[[idc]]))
go_col  <- as.character(pg[[goc]])
for (i in seq_along(pid_col)) {
  if (is.na(go_col[i]) || !nzchar(trimws(go_col[i]))) next
  ids <- unique(regmatches(go_col[i], gregexpr("GO:\\d{7}", go_col[i]))[[1]])
  if (length(ids) == 0) next
  pid <- sub("\\..*$", "", pid_col[i])                 # 去版本号(与 crosswalk 键一致)
  g <- unname(pwa2gene[pid])
  if (length(g) != 1 || is.na(g) || !(g %in% anno$gene_id)) next   # 只保留本研究基因
  ## go_term/namespace 先留空: eggNOG/InterProScan 的 Description 是蛋白描述而非 GO 术语名,
  ## 故术语名统一由下方 GO_OBO(go-basic.obo)补全; 无 OBO 时回退为 GO id。
  rows[[length(rows)+1]] <- data.frame(gene_id = g, go_id = ids,
                                       go_term = rep(NA_character_, length(ids)),
                                       namespace = NA_character_, stringsAsFactors = FALSE)
}
g <- if (length(rows)) unique(do.call(rbind, rows)) else
     data.frame(gene_id=character(0), go_id=character(0), go_term=character(0),
                namespace=character(0), stringsAsFactors=FALSE)

## 可选: 用 go-basic.obo 补 GO term 名称与命名空间
obo <- Sys.getenv("GO_OBO", "~/eggnog_data/go-basic.obo")  # 默认用持久 obo 补术语名/命名空间
if (nzchar(obo) && file.exists(obo) && nrow(g) > 0) {
  cat(sprintf("[annotation] 用 GO OBO 补术语/命名空间: %s\n", obo))
  L <- readLines(obo, warn = FALSE)
  id <- grep("^id: GO:", L); nm <- grep("^name:", L); ns <- grep("^namespace:", L)
  gid <- sub("^id: ", "", L[id]); gnm <- sub("^name: ", "", L[nm]); gns <- sub("^namespace: ", "", L[ns])
  ## 每个 [Term] 块的顺序为 id,name,...,namespace -> 按出现顺序对齐(取每个 id 后最近的 name/ns)
  name_map <- setNames(gnm, sub(".*", "", gnm)); name_map <- gnm; names(name_map) <- gid
  ns_map <- gns; names(ns_map) <- gid
  m <- match(g$go_id, gid)
  g$go_term <- ifelse(is.na(g$go_term) | !nzchar(g$go_term), ifelse(is.na(m), g$go_id, gnm[pmax(m,1)]), g$go_term)
  nsv <- ifelse(is.na(m), "NA", gns[pmax(m,1)])
  g$namespace <- ifelse(grepl("biological_process", nsv), "BP",
                 ifelse(grepl("cellular_component", nsv), "CC",
                 ifelse(grepl("molecular_function", nsv), "MF", "NA")))
}
miss <- is.na(g$go_term) | !nzchar(trimws(g$go_term)); g$go_term[miss] <- g$go_id[miss]
g$namespace[is.na(g$namespace) | !nzchar(g$namespace)] <- "NA"

write.table(g, out_go, sep = "\t", quote = FALSE, row.names = FALSE,
            col.names = c("gene_id","go_id","go_term","namespace"))
cat(sprintf("[annotation] 生成 GO 文件: %s\n  基因 %d 个, gene-term %d 条, 唯一 GO %d\n",
            out_go, length(unique(g$gene_id)), nrow(g), length(unique(g$go_id))))
## ---- 3) 基因->KO(供 KEGG 通路富集): 复用 eggNOG KEGG_ko 列 + PROT->gene 桥接 ----
##   KEGG 通路定义(KO->pathway、pathway 名)另由 KEGG 文件提供(ko2pathway.tsv / pathway_names.tsv)
koc <- NA_integer_
## 优先按表头名 KEGG_ko 定位(eggNOG); 否则按 "ko:K#####" 内容定位(KEGG_ko 列特征, 避免误配 seed_ortholog)
hdr_line <- grep("^#query", raw_lines, value = TRUE)
if (length(hdr_line)) { coln <- strsplit(sub("^#", "", hdr_line[1]), "\t", fixed = TRUE)[[1]]
  koc <- grep("^KEGG_ko$", coln)[1] }
if (is.na(koc) || koc > ncol(pg)) {
  for (j in seq_len(ncol(pg))) { v <- head(as.character(pg[[j]]), 800)
    if (any(grepl("ko:K[0-9]{5}", v[!is.na(v)]))) { koc <- j; break } }
}
ko_out <- Sys.getenv("KO_OUT", "~/eggnog_data/ko_by_gene.tsv")
if (!is.na(koc)) {
  ko_col <- as.character(pg[[koc]])
  krows <- list()
  for (i in seq_along(pid_col)) {
    if (is.na(ko_col[i]) || !nzchar(trimws(ko_col[i]))) next
    kos <- unique(regmatches(ko_col[i], gregexpr("K[0-9]{5}", ko_col[i]))[[1]])
    if (length(kos) == 0) next
    pid <- sub("\\..*$", "", pid_col[i])
    gg <- unname(pwa2gene[pid])
    if (length(gg) != 1 || is.na(gg) || !(gg %in% anno$gene_id)) next
    krows[[length(krows)+1]] <- data.frame(gene_id = gg, ko_id = kos, stringsAsFactors = FALSE)
  }
  kb <- if (length(krows)) unique(do.call(rbind, krows)) else
        data.frame(gene_id=character(0), ko_id=character(0), stringsAsFactors=FALSE)
  write.table(kb, ko_out, sep="\t", quote=FALSE, row.names=FALSE, col.names=c("gene_id","ko_id"))
  cat(sprintf("[annotation] 生成 基因->KO: %s (基因 %d, gene-KO %d, KO %d)\n",
              ko_out, length(unique(kb$gene_id)), nrow(kb), length(unique(kb$ko_id))))
} else cat("[annotation] 警告: 未找到 KEGG_ko 列, 跳过基因->KO(KEGG 富集将不可用)\n")

cat(sprintf("\n[下一步] GO_ANNOT_FILE=%s Rscript Pipeline/06_go.R\n", out_go))
cat("         然后 Rscript Pipeline/07_enrichment.R (富集依赖 06_go.rds)\n")
cat("KEGG 富集需 ~/eggnog_data/{ko_by_gene.tsv, ko2pathway.tsv, pathway_names.tsv} 就位; 就位则 07 自动做 KEGG。\n")
