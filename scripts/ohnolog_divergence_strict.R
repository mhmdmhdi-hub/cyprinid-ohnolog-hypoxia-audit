library(ggplot2)

df <- read.csv("~/carp_hypoxia/deseq2_annotated.csv", stringsAsFactors = FALSE)

# 1. Keep annotated genes on A/B chromosomes with valid padj
df <- df[!is.na(df$padj) & !is.na(df$gene_name) & df$gene_name != "NA", ]
df <- df[grepl("^[AB][0-9]+$", df$chr), ]

df$sub <- substr(df$chr, 1, 1)
df$chr_num <- sub("^[AB]", "", df$chr)

# 2. Keep strict 1:1 A-B pairs on matching homeologous chromosomes (A_i vs B_i)
a_df <- df[df$sub == "A", ]
b_df <- df[df$sub == "B", ]

a_counts <- table(a_df$gene_name)
b_counts <- table(b_df$gene_name)

valid_names <- intersect(names(a_counts[a_counts == 1]), names(b_counts[b_counts == 1]))
a_1to1 <- a_df[a_df$gene_name %in% valid_names, c("gene_name", "gene_id", "chr", "chr_num", "log2FoldChange", "padj")]
b_1to1 <- b_df[b_df$gene_name %in% valid_names, c("gene_name", "gene_id", "chr", "chr_num", "log2FoldChange", "padj")]

pairs <- merge(a_1to1, b_1to1, by = "gene_name", suffixes = c("_A", "_B"))
strict <- pairs[pairs$chr_num_A == pairs$chr_num_B, ]

# 3. Compute divergence metrics
strict$delta_lfc <- strict$log2FoldChange_A - strict$log2FoldChange_B
strict$abs_delta <- abs(strict$delta_lfc)

a_resp <- strict$padj_A < 0.05 & abs(strict$log2FoldChange_A) >= 1
b_resp <- strict$padj_B < 0.05 & abs(strict$log2FoldChange_B) >= 1
any_resp <- a_resp | b_resp

strict$category <- "Non-responsive"
strict$category[any_resp & strict$abs_delta < 1.5] <- "Conserved"
strict$category[any_resp & strict$abs_delta >= 1.5] <- "Diverged (|dLFC| >= 1.5)"
strict$category[a_resp & b_resp & (strict$log2FoldChange_A * strict$log2FoldChange_B < 0)] <- "Opposite"

# Save strict table
write.csv(strict, "~/carp_hypoxia/strict_homeolog_pairs.csv", row.names = FALSE)

# 4. Print summary stats
cat("\n=== Strict Homeologous Pairs (A_i vs B_i) ===\n")
cat("Total 1:1 A-B pairs (all chromosomes):", nrow(pairs), "\n")
cat("Strict matching A_i/B_i pairs:        ", nrow(strict), "\n\n")
print(table(strict$category))

cor_all  <- cor(strict$log2FoldChange_A, strict$log2FoldChange_B, method = "pearson")
cor_resp <- cor(strict$log2FoldChange_A[any_resp], strict$log2FoldChange_B[any_resp], method = "pearson")
cat(sprintf("\nPearson r (all strict pairs):        %.3f\n", cor_all))
cat(sprintf("Pearson r (responsive pairs only):   %.3f\n\n", cor_resp))

cat("=== Top 10 Strictly Diverged Homeologs (by |delta LFC|) ===\n")
top_div <- strict[strict$category %in% c("Diverged (|dLFC| >= 1.5)", "Opposite"), ]
top_div <- top_div[order(-top_div$abs_delta), ]
print(head(top_div[, c("gene_name", "chr_A", "log2FoldChange_A", "chr_B", "log2FoldChange_B", "delta_lfc", "category")], 10), row.names = FALSE)

# 5. Scatter plot
p <- ggplot(strict, aes(x = log2FoldChange_A, y = log2FoldChange_B, color = category)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "black") +
  geom_point(alpha = 0.6, size = 1.6) +
  scale_color_manual(values = c(
    "Non-responsive" = "grey80",
    "Conserved" = "#2b83ba",
    "Diverged (|dLFC| >= 1.5)" = "#fdae61",
    "Opposite" = "#d7191c"
  )) +
  theme_bw(base_size = 13) +
  labs(
    title = "Ohnolog Expression Divergence in C. carpio Brain (Hypoxia)",
    subtitle = sprintf("Strict A_i/B_i pairs (n = %d) | Responsive Pearson r = %.2f", nrow(strict), cor_resp),
    x = "Subgenome A log2FoldChange",
    y = "Subgenome B log2FoldChange",
    color = "Pattern"
  )

ggsave("~/carp_hypoxia/ohnolog_scatter_strict.png", plot = p, width = 8, height = 6, dpi = 300)
cat("\nScatter plot saved to ~/carp_hypoxia/ohnolog_scatter_strict.png\n")
