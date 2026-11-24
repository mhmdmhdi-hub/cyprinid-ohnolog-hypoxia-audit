#!/usr/bin/env Rscript
workdir   <- path.expand("~/carassius_native_2025")
carp_file <- path.expand("~/carp_hypoxia/deseq2_carp_brain_GENE_level.csv")
q         <- read.csv(file.path(workdir, "native_verified_quartets.csv"), stringsAsFactors = FALSE)
carp      <- read.csv(carp_file, stringsAsFactors = FALSE)
ccar      <- read.csv(file.path(workdir, "carassius_native_deseq2.csv"), stringsAsFactors = FALSE)
ccar_sex  <- read.csv(file.path(workdir, "carassius_native_deseq2_sex_adj.csv"), stringsAsFactors = FALSE)
norm_ccar <- read.csv(file.path(workdir, "carassius_native_norm_counts.csv"), stringsAsFactors = FALSE)
rownames(norm_ccar) <- norm_ccar$gene_id

merge_quartets <- function(cc_df) {
  m <- merge(q, carp[, c("gene_id","baseMean","log2FoldChange","lfcSE","padj")], by.x = "carp_gene_A", by.y = "gene_id")
  names(m)[(ncol(m)-3):ncol(m)] <- c("carp_bm_A","carp_lfc_A","carp_se_A","carp_padj_A")
  m <- merge(m, carp[, c("gene_id","baseMean","log2FoldChange","lfcSE","padj")], by.x = "carp_gene_B", by.y = "gene_id")
  names(m)[(ncol(m)-3):ncol(m)] <- c("carp_bm_B","carp_lfc_B","carp_se_B","carp_padj_B")
  m <- merge(m, cc_df[, c("gene_id","baseMean","log2FoldChange","lfcSE","padj")], by.x = "ccar_orth_of_carpA", by.y = "gene_id")
  names(m)[(ncol(m)-3):ncol(m)] <- c("ccar_bm_A","ccar_lfc_A","ccar_se_A","ccar_padj_A")
  m <- merge(m, cc_df[, c("gene_id","baseMean","log2FoldChange","lfcSE","padj")], by.x = "ccar_orth_of_carpB", by.y = "gene_id")
  names(m)[(ncol(m)-3):ncol(m)] <- c("ccar_bm_B","ccar_lfc_B","ccar_se_B","ccar_padj_B")
  m <- m[!is.na(m$carp_padj_A) & !is.na(m$carp_padj_B) & !is.na(m$ccar_padj_A) & !is.na(m$ccar_padj_B), ]
  m$carp_dlfc <- m$carp_lfc_A - m$carp_lfc_B
  m$ccar_dlfc <- m$ccar_lfc_A - m$ccar_lfc_B
  m$carp_resp <- (m$carp_padj_A < 0.05 & abs(m$carp_lfc_A) >= 1) | (m$carp_padj_B < 0.05 & abs(m$carp_lfc_B) >= 1)
  m$ccar_resp <- (m$ccar_padj_A < 0.05 & abs(m$ccar_lfc_A) >= 1) | (m$ccar_padj_B < 0.05 & abs(m$ccar_lfc_B) >= 1)
  m$carp_opp  <- (m$carp_padj_A < 0.05 & abs(m$carp_lfc_A) >= 1) & (m$carp_padj_B < 0.05 & abs(m$carp_lfc_B) >= 1) & (m$carp_lfc_A * m$carp_lfc_B < 0)
  m$ccar_opp  <- (m$ccar_padj_A < 0.05 & abs(m$ccar_lfc_A) >= 1) & (m$ccar_padj_B < 0.05 & abs(m$ccar_lfc_B) >= 1) & (m$ccar_lfc_A * m$ccar_lfc_B < 0)
  m$carp_div  <- m$carp_resp & (abs(m$carp_dlfc) >= 1.5 | m$carp_opp)
  m$ccar_div  <- m$ccar_resp & (abs(m$ccar_dlfc) >= 1.5 | m$ccar_opp)
  m
}

run_fisher <- function(df, label, filt) {
  sub <- df[filt, ]; co <- sub[sub$carp_resp & sub$ccar_resp, ]
  sd <- sum(co$carp_div & co$ccar_div); cp <- sum(co$carp_div & !co$ccar_div)
  cc <- sum(!co$carp_div & co$ccar_div); bc <- sum(!co$carp_div & !co$ccar_div)
  ft <- fisher.test(matrix(c(sd, cc, cp, bc), nrow = 2), alternative = "two.sided")
  data.frame(analysis = label, background = nrow(sub), co_responsive = nrow(co),
             shared_div = sd, carp_only = cp, ccar_only = cc, both_cons = bc,
             odds_ratio = unname(ft$estimate), ci_low = ft$conf.int[1], ci_high = ft$conf.int[2],
             p_value = ft$p.value, stringsAsFactors = FALSE)
}

m_main <- merge_quartets(ccar); m_sex <- merge_quartets(ccar_sex)
summary_tbl <- rbind(
  run_fisher(m_main, "Primary (~condition, max baseMean >= 50)", pmax(m_main$carp_bm_A, m_main$carp_bm_B) >= 50 & pmax(m_main$ccar_bm_A, m_main$ccar_bm_B) >= 50),
  run_fisher(m_sex,  "Sex-adjusted (~sex+condition, max baseMean >= 50)", pmax(m_sex$carp_bm_A, m_sex$carp_bm_B) >= 50 & pmax(m_sex$ccar_bm_A, m_sex$ccar_bm_B) >= 50),
  run_fisher(m_main, "Strict both-copy (min baseMean >= 10)", pmin(m_main$carp_bm_A, m_main$carp_bm_B) >= 10 & pmin(m_main$ccar_bm_A, m_main$ccar_bm_B) >= 10),
  run_fisher(m_main, "Strict both-copy (min baseMean >= 20)", pmin(m_main$carp_bm_A, m_main$carp_bm_B) >= 20 & pmin(m_main$ccar_bm_A, m_main$ccar_bm_B) >= 20)
)
write.csv(summary_tbl, "results/sensitivity_fisher_summary.csv", row.names = FALSE)

co_main <- m_main[pmax(m_main$carp_bm_A, m_main$carp_bm_B) >= 50 & pmax(m_main$ccar_bm_A, m_main$ccar_bm_B) >= 50 & m_main$carp_resp & m_main$ccar_resp, ]
samp_cols <- c("N14","N15","N16","N17","A1","A2","A3","A4")
cond_vec  <- factor(c(rep("Normoxia",4), rep("Anoxia",4)), levels = c("Normoxia","Anoxia"))
sex_vec   <- factor(c("F","M","F","M","M","M","M","F"))
p_mat <- t(mapply(function(gA, gB) {
  ratio <- log2(as.numeric(norm_ccar[gA, samp_cols]) + 1) - log2(as.numeric(norm_ccar[gB, samp_cols]) + 1)
  c(p_int = summary(lm(ratio ~ cond_vec))$coefficients["cond_vecAnoxia", "Pr(>|t|)"],
    p_int_sex = summary(lm(ratio ~ sex_vec + cond_vec))$coefficients["cond_vecAnoxia", "Pr(>|t|)"])
}, co_main$ccar_orth_of_carpA, co_main$ccar_orth_of_carpB))
co_main$ccar_p_paired_int     <- p_mat[, "p_int"]
co_main$ccar_fdr_paired_int   <- p.adjust(co_main$ccar_p_paired_int, method = "BH")
co_main$ccar_p_paired_int_sex <- p_mat[, "p_int_sex"]
write.csv(co_main, "results/native_coresponsive_full.csv", row.names = FALSE)
