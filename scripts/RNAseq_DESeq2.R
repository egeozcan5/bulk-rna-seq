# Run from the project folder containing data/ and scripts/.
library(DESeq2)
library(org.Hs.eg.db)
library(AnnotationDbi)

output_dir <- "results/DESeq2"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Counts and sample groups
counts_raw <- read.delim("data/gene_counts.txt", comment.char = "#")
count_matrix <- counts_raw[, 7:ncol(counts_raw)]
rownames(count_matrix) <- counts_raw$Geneid

sample_info <- data.frame(
  row.names = colnames(count_matrix),
  condition = factor(
    c("Proliferating", "Proliferating", "Proliferating",
      "Early", "Early", "Early", "Deep", "Deep", "Deep"),
    levels = c("Proliferating", "Early", "Deep")
  )
)

# Fit the model.
dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = sample_info,
  design = ~ condition
)
dds <- DESeq(dds)

# PCA using variance-stabilized counts.
vsd <- vst(dds, blind = TRUE)
pdf(file.path(output_dir, "PCA.pdf"), width = 8, height = 6)
print(plotPCA(vsd, intgroup = "condition"))
dev.off()

# Early vs Proliferating.
res_early <- results(dds, contrast = c("condition", "Early", "Proliferating"))
res_early <- res_early[order(res_early$padj), ]
pdf(file.path(output_dir, "Early_vs_Proliferating_MA.pdf"), width = 8, height = 6)
plotMA(res_early, ylim = c(-10, 10))
dev.off()

sig_early <- subset(res_early, padj < 0.05 & abs(log2FoldChange) > 1)
gene_ids_early <- sub("\\..*$", "", rownames(sig_early))
symbols_early <- mapIds(
  org.Hs.eg.db,
  keys = gene_ids_early,
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)

report_early <- data.frame(
  Gene = unname(symbols_early),
  Ensembl_ID = rownames(sig_early),
  as.data.frame(sig_early),
  row.names = NULL
)
report_early <- report_early[order(report_early$padj, -abs(report_early$log2FoldChange)), ]
up_early <- subset(report_early, log2FoldChange > 1)
down_early <- subset(report_early, log2FoldChange < -1)

write.csv(report_early, file.path(output_dir, "Early_vs_Proliferating_DEGs.csv"), row.names = FALSE)
write.csv(up_early, file.path(output_dir, "Early_vs_Proliferating_Upregulated.csv"), row.names = FALSE)
write.csv(down_early, file.path(output_dir, "Early_vs_Proliferating_Downregulated.csv"), row.names = FALSE)

report_early$Regulation <- ifelse(report_early$log2FoldChange > 0, "Up", "Down")
report_early$Comparison <- rep("Early_vs_Proliferating", nrow(report_early))

# Deep vs Proliferating.
res_deep <- results(dds, contrast = c("condition", "Deep", "Proliferating"))
pdf(file.path(output_dir, "Deep_vs_Proliferating_MA.pdf"), width = 8, height = 6)
plotMA(res_deep)
dev.off()

sig_deep <- subset(res_deep, padj < 0.05 & abs(log2FoldChange) > 1)
gene_ids_deep <- sub("\\..*$", "", rownames(sig_deep))
symbols_deep <- mapIds(
  org.Hs.eg.db,
  keys = gene_ids_deep,
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)

report_deep <- data.frame(
  Gene = unname(symbols_deep),
  Ensembl_ID = rownames(sig_deep),
  as.data.frame(sig_deep),
  row.names = NULL
)
report_deep <- report_deep[order(report_deep$padj, -abs(report_deep$log2FoldChange)), ]
up_deep <- subset(report_deep, log2FoldChange > 1)
down_deep <- subset(report_deep, log2FoldChange < -1)

write.csv(report_deep, file.path(output_dir, "Deep_vs_Proliferating_DEGs.csv"), row.names = FALSE)
write.csv(up_deep, file.path(output_dir, "Deep_vs_Proliferating_Upregulated.csv"), row.names = FALSE)
write.csv(down_deep, file.path(output_dir, "Deep_vs_Proliferating_Downregulated.csv"), row.names = FALSE)

report_deep$Regulation <- ifelse(report_deep$log2FoldChange > 0, "Up", "Down")
report_deep$Comparison <- rep("Deep_vs_Proliferating", nrow(report_deep))

# Combine the two comparisons; a gene can appear in both.
all_DEGs <- rbind(report_early, report_deep)
all_DEGs <- all_DEGs[order(all_DEGs$padj), c(
  "Comparison", "Regulation", "Gene", "Ensembl_ID",
  "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"
)]
write.csv(all_DEGs, file.path(output_dir, "All_Differentially_Expressed_Genes.csv"),
          row.names = FALSE)
