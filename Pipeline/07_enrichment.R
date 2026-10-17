# ============================================================
# 07_enrichment.R  —  PFAM 结构域 + TF 家族 + GO 超几何富集
#   GO 映射由前一步 06_go.rds 提供(不再在此抓取); 本步只做富集检验。
#   超几何(PFAM/TF/GO/KEGG) vs **"有资格被判为 DEG"的背景**(逐分析取 padj 非 NA 的基因),
#   BH 校正; 特异性用 2x2 Fisher; 另用旧背景(全部过滤基因)做敏感性对照。
#   方法与文献见 tables/enrichment_methods.txt, 对照结果见 tables/enrich_background_sensitivity.csv。
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
count_mat <- f$count_mat
d  <- readRDS(file.path(rdsdir, "04_deseq.rds"))
tg <- readRDS(file.path(rdsdir, "05_targets.rds"))

## 7.1 从描述列提取 PFAM 域
extract_pfam <- function(desc) {
  if (is.na(desc) || desc == "") return(character(0))
  parts <- strsplit(desc, "&&")[[1]]
  pfs <- character(0)
  for (p in parts) { m <- regmatches(p, gregexpr("PF[0-9]{5}", p))[[1]]; if (length(m)) pfs <- c(pfs, m) }
  unique(pfs)
}
if (!"gene_description" %in% colnames(anno))
  stop("注释表缺少 'gene_description' 列, 无法提取 PFAM 域; 请在输入中提供该列。")
pfam_list <- setNames(lapply(anno$gene_description, extract_pfam), anno$gene_id)

## 7.2 TF 家族列表 (lapply 逐元素构造, 避免 ifelse 对 list 的已知缺陷)
fam_list <- setNames(lapply(seq_along(anno$Family), function(i)
  if (anno$Family[i] == "-") character(0) else anno$Family[i]), anno$gene_id)

## 7.3 通用超几何富集 (向量化 + phyper 参数校验)
enrich_hyper <- function(test_ids, bg_ids, feature_list, feat_name, min_bg = 5) {
  test_ids <- intersect(test_ids, names(feature_list))
  bg_ids   <- intersect(bg_ids,   names(feature_list))
  K <- length(test_ids); N <- length(bg_ids)
  if (K == 0 || N == 0) return(data.frame())
  test_feats <- unlist(feature_list[test_ids], use.names = FALSE)
  bg_feats   <- unlist(feature_list[bg_ids],   use.names = FALSE)
  if (length(test_feats) == 0) return(data.frame())
  tt <- table(test_feats); bb <- table(bg_feats)
  allfeat <- union(names(tt), names(bb))
  k <- as.integer(tt[allfeat]); k[is.na(k)] <- 0
  m <- as.integer(bb[allfeat]); m[is.na(m)] <- 0
  keep <- k > 0 & m >= min_bg
  k <- k[keep]; m <- m[keep]; allfeat <- allfeat[keep]
  if (length(k) == 0) return(data.frame())
  valid <- (k - 1 >= 0) & (m >= 0) & (N - m >= 0) & (K >= 0) & (k - 1 <= K)
  p <- rep(1, length(k))
  p[valid] <- phyper(k[valid] - 1, m[valid], N - m[valid], K, lower.tail = FALSE)
  res <- data.frame(feature = allfeat, in_test = k, in_bg = m,
                    test_total = K, bg_total = N, pvalue = p, stringsAsFactors = FALSE)
  res$padj <- p.adjust(res$pvalue, method = "BH")
  res$type <- feat_name
  res[order(res$pvalue), ]
}

## 7.3b 背景(检验总体)定义 —— 必须是"有机会被判为 DEG 的基因"
##   超几何检验把 DEG 集看作从总体中抽出的样本, 要求总体内每个基因都有非零的入选机会。
##   DESeq2 的独立过滤/Cook's 距离会把部分基因 padj 置 NA, 这些基因不可能进入 DEG 集
##   (这些基因通常表达量更低), 把它们留在背景里等于
##   让"入选概率恒为 0"的基因参与抽样, 破坏可交换性; 偏倚方向还随"表达-功能"相关性变化。
##   故逐分析使用各自的合格背景:
##     单对比          -> 该对比 padj 非 NA 的基因
##     DEG 并集        -> 三个对比合格基因的并集 (任一对比可判即可入选)
##     靶标集(双对比)  -> OE 与 RNAi 两个对比都合格的交集
##   旧的"全部过滤基因"口径保留为敏感性对照 (tables/enrich_background_sensitivity.csv)。
elig <- function(tab) tab$gene_id[!is.na(tab$padj)]
bg_oe     <- elig(d$t_OE_CK)
bg_ri     <- elig(d$t_RNAi_CK)
bg_or     <- elig(d$t_OE_RNAi)
bg_union  <- union(union(bg_oe, bg_ri), bg_or)
bg_target <- intersect(bg_oe, bg_ri)
bg_all    <- rownames(count_mat)          # 旧口径(仅作敏感性对照)
union_deg <- union(union(d$deg_OE, d$deg_RNAi), d$deg_OR)
cat(sprintf("[背景] 有资格判 DEG 的基因数: OE_vs_CK=%d, RNAi_vs_CK=%d, OE_vs_RNAi=%d, 并集=%d, 靶标(OE∩RNAi)=%d (全部过滤基因=%d)\n",
            length(bg_oe), length(bg_ri), length(bg_or), length(bg_union),
            length(bg_target), length(bg_all)))

## 包装: 每个检验算两次(合格背景=主结果; 全部过滤基因=敏感性对照), 便于做背景稳健性说明
sens <- list()
enrich_main <- function(key, ids, bg_ids, feature_list, feat_name) {
  main <- enrich_hyper(ids, bg_ids, feature_list, feat_name)
  alt  <- enrich_hyper(ids, bg_all,  feature_list, feat_name)
  smry <- function(r) {
    if (is.null(r) || nrow(r) == 0)
      return(list(n = 0L, nsig = 0L, minpadj = NA_real_, feat = NA_character_, K = 0L))
    i <- which.min(r$pvalue)
    list(n = nrow(r), nsig = sum(r$padj < 0.05, na.rm = TRUE),
         minpadj = min(r$padj, na.rm = TRUE), feat = as.character(r$feature[i]),
         K = r$test_total[1])
  }
  a <- smry(main); b <- smry(alt)
  sens[[key]] <<- data.frame(
    table = key, test_set_size = a$K, n_features = a$n,
    min_padj_tested_bg = a$minpadj,  n_sig_tested_bg = a$nsig,
    min_padj_allfiltered_bg = b$minpadj, n_sig_allfiltered_bg = b$nsig,
    top_feature = a$feat, stringsAsFactors = FALSE)
  if (nrow(main) > 0)
    main$background <- "genes eligible for DEG calling (padj non-NA)"
  main
}

enr <- list()
enr[["PFAM_DEG_union"]]     <- enrich_main("PFAM_DEG_union",     union_deg, bg_union,  pfam_list, "PFAM")
enr[["PFAM_activated"]]     <- enrich_main("PFAM_activated",     tg$act,    bg_target, pfam_list, "PFAM")
enr[["PFAM_repressed"]]     <- enrich_main("PFAM_repressed",     tg$rep,    bg_target, pfam_list, "PFAM")
enr[["TFfamily_DEG_union"]] <- enrich_main("TFfamily_DEG_union", union_deg, bg_union,  fam_list,  "TFfamily")
enr[["TFfamily_activated"]] <- enrich_main("TFfamily_activated", tg$act,    bg_target, fam_list,  "TFfamily")
enr[["TFfamily_repressed"]] <- enrich_main("TFfamily_repressed", tg$rep,    bg_target, fam_list,  "TFfamily")

## 7.4 GO 富集 (GO 映射来自 06_go.rds; BH 校正; 过滤稀有项)
go_rds <- file.path(rdsdir, "06_go.rds")
go_ann <- if (file.exists(go_rds)) readRDS(go_rds) else NULL
go_status <- if (is.null(go_ann)) {
  "missing"
} else if (is.null(go_ann$status)) {
  "ok"
} else {
  go_ann$status
}
go_ok <- identical(go_status, "ok") && !is.null(go_ann$go_term2gene) &&
         nrow(go_ann$go_term2gene) > 0
go_stats <- NULL
if (go_ok) {
  gtg <- go_ann$go_term2gene
  gtg <- gtg[!duplicated(gtg[, c("gene_id", "go_id")]), ]   # 同基因同 GO 不重复计数
  go_list <- split(gtg$go_id, gtg$gene_id)
  go_meta <- gtg[!duplicated(gtg$go_id), c("go_id", "go_term", "namespace")]
  go_meta$go_term <- trimws(gsub("^[;[:space:]]+", "", go_meta$go_term))
  go_stats <- list(source = go_ann$source, source_file = go_ann$source_file, status = go_status,
                   genes = length(unique(gtg$gene_id)), terms = length(unique(gtg$go_id)),
                   gene_term = nrow(gtg),
                   terms_no_ns = sum(!go_meta$namespace %in% c("BP", "CC", "MF")))
  enrich_go <- function(key, ids, bg_ids) {
    r <- enrich_main(key, ids, bg_ids, go_list, "GO")
    if (nrow(r) == 0) return(r)
    r <- merge(r, go_meta, by.x = "feature", by.y = "go_id", all.x = TRUE)
    r[order(r$pvalue), ]
  }
  enr[["GO_DEG_union"]] <- enrich_go("GO_DEG_union", union_deg,  bg_union)
  enr[["GO_OE"]]        <- enrich_go("GO_OE",        d$deg_OE,   bg_oe)
  enr[["GO_RNAi"]]      <- enrich_go("GO_RNAi",      d$deg_RNAi, bg_ri)
  enr[["GO_OE_RNAi"]]   <- enrich_go("GO_OE_RNAi",   d$deg_OR,   bg_or)
  enr[["GO_activated"]] <- enrich_go("GO_activated", tg$act,     bg_target)
  enr[["GO_repressed"]] <- enrich_go("GO_repressed", tg$rep,     bg_target)
  cat(sprintf("[GO] 显著项(padj<0.05): DEG_union=%d, activated=%d, repressed=%d\n",
              sum(enr[["GO_DEG_union"]]$padj < 0.05, na.rm = TRUE),
              sum(enr[["GO_activated"]]$padj < 0.05,  na.rm = TRUE),
              sum(enr[["GO_repressed"]]$padj < 0.05,  na.rm = TRUE)))
} else {
  ## 必须显式告警: 空注释 -> 本步不产出 GO 富集, 但这**不等于**"GO 无显著富集"
  warning("[GO] 06 步骤未产出可用注释 (status=", go_status,
          "): 本次结果不含 GO 富集 —— 不要解读为'无富集'。请先修复 06_go.R 后重跑 06/07。")
}

## 7.5 KEGG 通路富集(目标物种原生 KO + KEGG 通路定义; 超几何 + BH, 同 GO/PFAM/TF 框架)
##   基因->KO: eggNOG 桥接(ko_by_gene.tsv); KO->pathway 与名称: KEGG 文件(ko2pathway/pathway_names)
kegg_dir  <- Sys.getenv("KEGG_DIR", "~/eggnog_data")
ko_gene_f <- file.path(kegg_dir, "ko_by_gene.tsv")
ko2path_f <- file.path(kegg_dir, "ko2pathway.tsv")
pname_f   <- file.path(kegg_dir, "pathway_names.tsv")
if (all(file.exists(ko_gene_f), file.exists(ko2path_f), file.exists(pname_f))) {
  kb <- read.delim(ko_gene_f, stringsAsFactors = FALSE)
  k2p <- read.delim(ko2path_f, header = FALSE, stringsAsFactors = FALSE, col.names = c("ko","pathway"))
  k2p$ko <- sub("^ko:", "", k2p$ko); k2p$pathway <- sub("^path:", "", k2p$pathway)
  k2p$pathway <- sub("^ko", "map", k2p$pathway)   # koXXXXX 与 mapXXXXX 是同一通路, 统一为 map 以去重并对齐名称表
  k2p <- unique(k2p)
  pname <- read.delim(pname_f, header = FALSE, stringsAsFactors = FALSE, col.names = c("pathway","name"))
  ko2paths  <- split(k2p$pathway, k2p$ko)                 # KO -> pathway(s)
  gene2path <- lapply(split(kb$ko_id, kb$gene_id), function(kos)
    unique(unlist(ko2paths[kos], use.names = FALSE)))    # 基因 -> pathway(s)
  gene2path <- gene2path[lengths(gene2path) > 0]
  kegg_stats <- list(genes = length(unique(kb$gene_id)), kos = length(unique(kb$ko_id)),
                     genes_mapped = length(gene2path),
                     pathways = length(unique(unlist(gene2path))))
  cat(sprintf("[KEGG] 基因->通路: %d 基因, %d 通路\n",
              length(gene2path), length(unique(unlist(gene2path)))))
  ## 外部数据库文件留痕 (路径/大小/md5/mtime): KEGG 定义文件不随仓库版本化
  kegg_prov <- c(sprintf("KEGG source files (%s)", Sys.time()),
                 unlist(lapply(c(ko_gene_f, ko2path_f, pname_f), function(p)
                   c(sprintf("  %s", p),
                     sprintf("    size=%s bytes; md5=%s; mtime=%s",
                             format(file.size(p), big.mark = ","),
                             unname(tools::md5sum(p)), format(file.mtime(p)))))))
  writeLines(c(kegg_prov,
               sprintf("  基因->KO: %d 基因 / %d KO; KO->通路表 %d 行; 通路名 %d 行",
                       length(unique(kb$gene_id)), length(unique(kb$ko_id)), nrow(k2p), nrow(pname)),
               "  提醒: 上述文件位于工程目录之外, 仅凭本仓库无法复现 KEGG 结果; 建议随补充材料归档"),
             file.path(outdir, "tables", "KEGG_source_provenance.txt"))

  enrich_kegg <- function(key, ids, bg_ids) {
    r <- enrich_main(key, ids, bg_ids, gene2path, "KEGG")
    if (nrow(r) == 0) return(r)
    r <- merge(r, pname, by.x = "feature", by.y = "pathway", all.x = TRUE)
    r[order(r$pvalue), ]
  }
  enr[["KEGG_DEG_union"]] <- enrich_kegg("KEGG_DEG_union", union_deg,  bg_union)
  enr[["KEGG_OE"]]        <- enrich_kegg("KEGG_OE",        d$deg_OE,   bg_oe)
  enr[["KEGG_RNAi"]]      <- enrich_kegg("KEGG_RNAi",      d$deg_RNAi, bg_ri)
  enr[["KEGG_OE_RNAi"]]   <- enrich_kegg("KEGG_OE_RNAi",   d$deg_OR,   bg_or)
  cat(sprintf("[KEGG] 显著(padj<0.05): DEG_union=%d OE=%d RNAi=%d OE_RNAi=%d\n",
              sum(enr[["KEGG_DEG_union"]]$padj<0.05,na.rm=TRUE),
              sum(enr[["KEGG_OE"]]$padj<0.05,na.rm=TRUE),
              sum(enr[["KEGG_RNAi"]]$padj<0.05,na.rm=TRUE),
              sum(enr[["KEGG_OE_RNAi"]]$padj<0.05,na.rm=TRUE)))
} else {
  kegg_stats <- NULL
  cat("[KEGG] 警告: 缺少 ko_by_gene.tsv / ko2pathway.tsv / pathway_names.tsv, 跳过 KEGG 富集\n")
}

## 7.6 特异性检验: 激活 vs 抑制 (2x2 Fisher), 过滤稀有域(min_total)
enrich_specific <- function(idsA, idsB, feature_list, feat_name, min_total = 3) {
  A <- intersect(idsA, names(feature_list)); B <- intersect(idsB, names(feature_list))
  nA <- length(A); nB <- length(B)
  if (nA == 0 || nB == 0) return(data.frame())
  ta <- table(unlist(feature_list[A], use.names = FALSE))
  tb <- table(unlist(feature_list[B], use.names = FALSE))
  allf <- union(names(ta), names(tb))
  a <- as.integer(ta[allf]); a[is.na(a)] <- 0
  b <- as.integer(tb[allf]); b[is.na(b)] <- 0
  keep <- (a + b) >= min_total
  a <- a[keep]; b <- b[keep]; allf <- allf[keep]
  if (length(allf) == 0) return(data.frame())
  c <- nA - a; dd <- nB - b
  p <- rep(1, length(a)); or <- rep(NA_real_, length(a))
  for (i in seq_along(a)) {
    ft <- fisher.test(matrix(c(a[i], c[i], b[i], dd[i]), nrow = 2), alternative = "two.sided")
    p[i] <- ft$p.value; or[i] <- ft$estimate
  }
  padj <- p.adjust(p, method = "BH")
  res <- data.frame(feature = allf, in_activated = a, in_repressed = b,
                    activated_total = nA, repressed_total = nB,
                    odds_ratio = or, pvalue = p, padj = padj, type = feat_name,
                    stringsAsFactors = FALSE)
  res$direction <- ifelse(a > b, "activated", ifelse(b > a, "repressed", "none"))
  res[order(res$pvalue), ]
}
enr[["PFAM_spec"]]     <- enrich_specific(tg$act, tg$rep, pfam_list, "PFAM")
enr[["TFfamily_spec"]] <- enrich_specific(tg$act, tg$rep, fam_list,  "TFfamily")

for (nm in names(enr))
  if (nrow(enr[[nm]]) > 0)
    write.csv(enr[[nm]], file.path(outdir, "tables", sprintf("enrich_%s.csv", nm)), row.names = FALSE)

## 7.7 背景口径敏感性对照 + 富集方法说明落盘
sens_df <- if (length(sens)) do.call(rbind, sens) else data.frame()
write.csv(sens_df, file.path(outdir, "tables", "enrich_background_sensitivity.csv"),
          row.names = FALSE)
n_changed <- if (nrow(sens_df)) sum(sens_df$n_sig_tested_bg != sens_df$n_sig_allfiltered_bg) else 0L

meth <- c(
  "Functional enrichment (TF RNA-seq) — step 07_enrichment.R",
  "============================================================",
  "Tests    : one-sided hypergeometric test (R phyper) for PFAM domains, TF families, GO terms and",
  "           KEGG pathways; Benjamini-Hochberg FDR (1995) applied within each feature type and gene set.",
  "           'Specificity' between activated and repressed target sets: 2x2 Fisher exact test.",
  "Feature sets: PFAM domains and TF families parsed from the input annotation table; GO from 06_go.rds;",
  "           KEGG via gene->KO (eggNOG) + KO->pathway (KEGG). Rare features (background count < 5) are",
  "           excluded because they cannot be tested reliably.",
  "",
  "Background (universe) — the genes that could have been called DE:",
  "  * DESeq2 independent filtering and Cook's-distance outlier removal set padj = NA for part of the",
  "    genes. Those genes cannot enter any DEG set and are typically lower-expressed, so keeping them",
  "    in the background inserts genes with a selection probability",
  "    of exactly zero, which violates the exchangeability assumption of the hypergeometric test.",
  "  * single contrast -> genes with non-NA padj in that contrast;",
  "  * DEG union       -> union of the three contrasts' testable genes;",
  "  * target sets     -> intersection of the OE-vs-CK and RNAi-vs-CK testable genes (a target must be",
  "                       callable in both contrasts);",
  "  * within each feature type the universe is further restricted to genes carrying that kind of",
  "    feature (reported as bg_total in every table).",
  sprintf("  universe sizes: OE_vs_CK=%d, RNAi_vs_CK=%d, OE_vs_RNAi=%d, union=%d, target=%d genes (filtered genes total=%d).",
          length(bg_oe), length(bg_ri), length(bg_or), length(bg_union), length(bg_target), length(bg_all)),
  "",
  "Sensitivity analysis: every test was repeated with the alternative background 'all filtered genes'",
  "  (the universe used before this revision). Both results are reported side by side in",
  "  tables/enrich_background_sensitivity.csv.",
  sprintf("  tables whose number of padj<0.05 features differs between the two backgrounds: %d of %d.",
          n_changed, nrow(sens_df)),
  "",
  "References (background choice):",
  "  Wijesooriya K, et al. Urgent need for consistent standards in functional enrichment analysis.",
  "    Brief Bioinform. 2022;23(3):bbac066. DOI:10.1093/bib/bbac066 | PMID:35272374",
  "  Timmons JA, Szkop KJ, Gallagher IJ. Multiple sources of bias confound functional enrichment",
  "    analysis of global -omics data. Genome Biol. 2015;16:186. DOI:10.1186/s13059-015-0761-7 | PMID:26346307",
  "  Benjamini Y, Hochberg Y. J R Stat Soc B. 1995;57(1):289-300.")
writeLines(meth, file.path(outdir, "tables", "enrichment_methods.txt"))

saveRDS(list(enr = enr, union_deg = union_deg,
             go_stats = go_stats, kegg_stats = kegg_stats,
             backgrounds = list(OE_vs_CK = bg_oe, RNAi_vs_CK = bg_ri, OE_vs_RNAi = bg_or,
                                union = bg_union, target = bg_target, all_filtered = bg_all),
             sens = sens_df),
        file.path(rdsdir, "07_enrichment.rds"))

cat(sprintf("[%s] 富集分析完成 (%d 个特征表; 背景口径敏感性: %d/%d 张表显著数不同, 见 tables/enrich_background_sensitivity.csv)\n",
            Sys.time(), length(enr), n_changed, nrow(sens_df)))
