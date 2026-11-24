#!/usr/bin/env Rscript
suppressPackageStartupMessages({ library(tximport); library(DESeq2) })
workdir <- path.expand("~/carassius_native_2025")
samples <- read.table("metadata/carassius_samples.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
samples$condition <- factor(samples$condition, levels = c("Normoxia", "Anoxia"))
samples$sex       <- factor(samples$sex, levels = c("female", "male"))
rownames(samples) <- samples$label
tx2gene <- read.table(file.path(workdir, "tx2gene_native.tsv"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
files   <- file.path(workdir, "salmon_quant", paste0(samples$srr, ".tabular")); names(files) <- samples$label
txi <- tximport(files, type = "salmon", tx2gene = tx2gene, ignoreTxVersion = FALSE, dropInfReps = TRUE)
dds <- DESeqDataSetFromTximport(txi, colData = samples, design = ~ condition)
keep <- rowSums(counts(dds) >= 10) >= 3
dds <- DESeq(dds[keep, ], quiet = TRUE)
res_unadj <- as.data.frame(results(dds, contrast = c("condition", "Anoxia", "Normoxia")))
res_unadj$gene_id <- rownames(res_unadj)
write.csv(res_unadj, file.path(workdir, "carassius_native_deseq2.csv"), row.names = FALSE)
norm_counts <- as.data.frame(counts(dds, normalized = TRUE)); norm_counts$gene_id <- rownames(norm_counts)
write.csv(norm_counts, file.path(workdir, "carassius_native_norm_counts.csv"), row.names = FALSE)
dds_sex <- DESeq(DESeqDataSetFromTximport(txi, colData = samples, design = ~ sex + condition)[keep, ], quiet = TRUE)
res_sex <- as.data.frame(results(dds_sex, contrast = c("condition", "Anoxia", "Normoxia")))
res_sex$gene_id <- rownames(res_sex)
write.csv(res_sex, file.path(workdir, "carassius_native_deseq2_sex_adj.csv"), row.names = FALSE)
