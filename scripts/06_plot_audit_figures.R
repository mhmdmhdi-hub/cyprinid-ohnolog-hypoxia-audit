#!/usr/bin/env Rscript
suppressPackageStartupMessages({ library(ggplot2); library(tximport); library(DESeq2) })
meta <- read.table("metadata/carassius_samples.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
meta$label <- factor(meta$label, levels = meta$label)

map_long <- rbind(
  data.frame(label = meta$label, Reference = "Non-native Goldfish (ASM336829v1)", MappingRate = meta$goldfish_map_pct),
  data.frame(label = meta$label, Reference = "Native C. carassius (2025, 50 chr)", MappingRate = meta$native_map_pct)
)
map_long$Reference <- factor(map_long$Reference, levels = c("Non-native Goldfish (ASM336829v1)", "Native C. carassius (2025, 50 chr)"))
p1 <- ggplot(map_long, aes(x = label, y = MappingRate, fill = Reference)) +
  geom_col(position = position_dodge(width = 0.72), width = 0.65, color = "grey20", linewidth = 0.25) +
  geom_text(aes(label = sprintf("%.1f%%", MappingRate)), position = position_dodge(width = 0.72), vjust = -0.45, size = 3.0) +
  scale_fill_manual(values = c("#bdbdbd", "#2b83ba")) +
  scale_y_continuous(limits = c(0, 76), breaks = seq(0, 70, 10)) +
  theme_bw(base_size = 12) + theme(legend.position = "top", panel.grid.minor = element_blank()) +
  labs(title = "Salmon Mapping Rate Gain on Native C. carassius Genome (2025)",
       subtitle = "PRJNA386629 brain libraries (n = 8) | Mean gain: +17.1% (46.98% -> 64.08%)",
       x = "Sample (N14-N17: Normoxia | A1-A4: Anoxia)", y = "Mapped Read Pairs (%)", fill = NULL)
ggsave("figures/01_mapping_rate_gain.png", p1, width = 8.5, height = 5.2, dpi = 300)

workdir <- path.expand("~/carassius_native_2025")
tx2gene <- read.table(file.path(workdir, "tx2gene_native.tsv"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
meta$condition <- factor(meta$condition, levels = c("Normoxia", "Anoxia"))
meta$sex <- factor(meta$sex, levels = c("female", "male"))
rownames(meta) <- meta$label
files <- file.path(workdir, "salmon_quant", paste0(meta$srr, ".tabular")); names(files) <- meta$label
txi <- tximport(files, type = "salmon", tx2gene = tx2gene, ignoreTxVersion = FALSE, dropInfReps = TRUE)
dds <- DESeqDataSetFromTximport(txi, colData = meta, design = ~ sex + condition)
dds <- dds[rowSums(counts(dds) >= 10) >= 3, ]
pca <- plotPCA(vst(dds, blind = TRUE), intgroup = c("condition", "sex"), returnData = TRUE)
pv  <- round(100 * attr(pca, "percentVar"), 1)
p2 <- ggplot(pca, aes(x = PC1, y = PC2, color = condition, shape = sex, label = name)) +
  geom_point(size = 4.2, stroke = 1.1) + geom_text(vjust = -1.1, size = 3.5, show.legend = FALSE) +
  scale_color_manual(values = c("Normoxia" = "#2b83ba", "Anoxia" = "#d7191c")) +
  scale_shape_manual(values = c("female" = 16, "male" = 17)) +
  theme_bw(base_size = 12) + theme(panel.grid.minor = element_blank()) +
  labs(title = "Native C. carassius Brain Transcriptome PCA (VST)",
       subtitle = "Condition dominates PC1 (89.3%); no sex-driven clustering (LFC r = 0.994)",
       x = sprintf("PC1: %.1f%% variance", pv[1]), y = sprintf("PC2: %.1f%% variance", pv[2]), color = "Condition", shape = "Sex")
ggsave("figures/02_native_vst_pca_sex.png", p2, width = 7.5, height = 5.2, dpi = 300)

sens <- read.csv("results/sensitivity_fisher_summary.csv", stringsAsFactors = FALSE)
sens$label <- factor(
  c("Primary (~condition)\nmax(baseMean) >= 50", "Sex-adjusted (~sex+condition)\nmax(baseMean) >= 50",
    "Strict both-copy\nmin(baseMean) >= 10", "Strict both-copy\nmin(baseMean) >= 20"),
  levels = rev(c("Primary (~condition)\nmax(baseMean) >= 50", "Sex-adjusted (~sex+condition)\nmax(baseMean) >= 50",
                 "Strict both-copy\nmin(baseMean) >= 10", "Strict both-copy\nmin(baseMean) >= 20")))
sens$sig <- ifelse(sens$p_value < 0.05, "p < 0.05 (fragile)", "Non-significant (null)")
sens$annot <- sprintf("OR = %.2f | p = %.4f | Shared = %d/%d", sens$odds_ratio, sens$p_value, sens$shared_div, sens$co_responsive)
p3 <- ggplot(sens, aes(x = odds_ratio, y = label, color = sig)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey45") +
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0.22, linewidth = 0.85) +
  geom_point(size = 3.8) + geom_text(aes(x = ci_high, label = annot), hjust = -0.08, size = 3.2, color = "grey20") +
  scale_x_log10(limits = c(0.12, 160), breaks = c(0.2, 0.5, 1, 2, 5, 10, 20)) +
  scale_color_manual(values = c("p < 0.05 (fragile)" = "#e66101", "Non-significant (null)" = "#5e3c99")) +
  theme_bw(base_size = 11.5) + theme(legend.position = "top", panel.grid.minor = element_blank()) +
  labs(title = "Sensitivity of Cross-Species Shared Ohnolog Divergence Enrichment",
       subtitle = "Enrichment collapses to null once both homeolog copies are required to be expressed",
       x = "Two-sided Fisher's Exact Odds Ratio (log10 scale, 95% CI)", y = NULL, color = NULL)
ggsave("figures/03_fisher_sensitivity_forest.png", p3, width = 9.8, height = 4.6, dpi = 300)
