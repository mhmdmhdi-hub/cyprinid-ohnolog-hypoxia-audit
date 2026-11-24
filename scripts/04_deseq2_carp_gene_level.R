library(tximport)
library(DESeq2)
library(ggplot2)

# ۱. تعریف مسیرها و جدول طراحی آزمایش
data_dir <- path.expand("~/carp_hypoxia/salmon_genes")
tx2gene  <- read.table("~/carp_hypoxia/tx2gene.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

samples <- data.frame(
  sample_id = c("SRR18578462", "SRR18578461", "SRR18578450",
                "SRR18578439", "SRR18578438", "SRR18578437"),
  label     = c("W1", "W2", "W3", "H1", "H2", "H3"),
  condition = factor(c("Control", "Control", "Control",
                       "Hypoxia", "Hypoxia", "Hypoxia"),
                     levels = c("Control", "Hypoxia"))
)
rownames(samples) <- samples$label

files <- file.path(data_dir, paste0(samples$sample_id, ".tabular"))
names(files) <- samples$label

# ۲. جمع‌بندی رونوشت‌ها به سطح ژن با tximport (با نادیده گرفتن پسوند .1)
txi <- tximport(files, type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE, dropInfReps = TRUE)

# ۳. ساخت آبجکت DESeq2 و فیلتر ژن‌های با بیان خیلی پایین
dds <- DESeqDataSetFromTximport(txi, colData = samples, design = ~ condition)
dds <- dds[rowSums(counts(dds)) >= 10, ]

# ۴. اجرای تحلیل آماری DESeq2
dds <- DESeq(dds)
res <- results(dds, contrast = c("condition", "Hypoxia", "Control"))
res_df <- as.data.frame(res)
res_df$gene_id <- rownames(res_df)
res_df <- res_df[order(res_df$padj), ]

# ذخیره جدول کامل نتایج در سطح ژن
write.csv(res_df, "~/carp_hypoxia/deseq2_carp_brain_GENE_level.csv", row.names = FALSE)

# ۵. چاپ خلاصه آماری در ترمینال
cat("\n=== خلاصه نتایج DESeq2 در سطح ژن (Hypoxia vs Control) ===\n")
cat("تعداد کل ژن‌های بررسی‌شده:", nrow(res_df), "\n")
cat("ژن‌های معنادار (padj < 0.05):", sum(res_df$padj < 0.05, na.rm = TRUE), "\n")
cat("  - افزایش بیان (Up, padj < 0.05 & log2FC > 1):",
    sum(res_df$padj < 0.05 & res_df$log2FoldChange > 1, na.rm = TRUE), "\n")
cat("  - کاهش بیان (Down, padj < 0.05 & log2FC < -1):",
    sum(res_df$padj < 0.05 & res_df$log2FoldChange < -1, na.rm = TRUE), "\n\n")

cat("=== ۱۰ ژن برتر (بر اساس padj) ===\n")
print(head(res_df[, c("gene_id", "baseMean", "log2FoldChange", "padj")], 10))

# ۶. رسم و ذخیره نمودار PCA در سطح ژن
vsd <- vst(dds, blind = FALSE)
pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))

p <- ggplot(pca_data, aes(PC1, PC2, color = condition, label = name)) +
  geom_point(size = 4) +
  geom_text(vjust = -1, size = 4, show.legend = FALSE) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_bw(base_size = 14) +
  ggtitle("PCA (Gene-level) - Cyprinus carpio Brain")

ggsave("~/carp_hypoxia/pca_brain_gene_level.png", plot = p, width = 7, height = 5, dpi = 300)
cat("\nنمودار PCA در مسیر ~/carp_hypoxia/pca_brain_gene_level.png ذخیره شد.\n")
