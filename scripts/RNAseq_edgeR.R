# Run from the repository root. Confirm sample order before running.
dir.create("results/edgeR", recursive = TRUE, showWarnings = FALSE)

# Load Required Packages
library(edgeR)
library(org.Hs.eg.db)
library(AnnotationDbi)

# Read featureCounts Output
counts_raw <- read.delim(
  "data/gene_counts.txt",
  comment.char = "#"
)

# Create Count Matrix
count_matrix <- counts_raw[, 7:ncol(counts_raw)]
rownames(count_matrix) <- counts_raw$Geneid

# Create Sample Information
sample_info <- data.frame(
  sample = colnames(count_matrix),
  condition = factor(c(
    "Proliferating",
    "Proliferating",
    "Proliferating",

    "Early",
    "Early",
    "Early",

    "Deep",
    "Deep",
    "Deep"

  ))
)

# Create DGEList Object
dge <- DGEList(counts = count_matrix, group = sample_info$condition)

# Filter Lowly Expressed Genes
keep <- filterByExpr(dge, group = sample_info$condition)

dge <- dge[keep, , keep.lib.sizes = FALSE]

# TMM Normalization
dge <- calcNormFactors(dge)

# QC - MDS Plot
plotMDS(
  dge,
  labels = sample_info$condition,
  col = c("blue","blue","blue", "orange","orange","orange", "red","red","red")
)

# Create Design Matrix
design <- model.matrix(~0 + sample_info$condition)
colnames(design) <- levels(sample_info$condition)

# Estimate Dispersion
dge <- estimateDisp(dge, design)
plotBCV(dge)

# Test for Differential Expression
fit <- glmQLFit(dge, design)
contrast_EP <- makeContrasts(Early - Proliferating, levels = design)

qlf_EP <- glmQLFTest(fit, contrast = contrast_EP)

# Gene Annotation
results_EP <- topTags(qlf_EP, n = Inf)$table

results_EP$Ensembl_ID <- rownames(results_EP)

gene_ids <- sub("\\..*$", "", results_EP$Ensembl_ID)

results_EP$Gene <- mapIds(
  org.Hs.eg.db,
  keys = gene_ids,
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)

results_EP <- results_EP[
  ,
  c("Gene", "Ensembl_ID", "logFC", "logCPM", "F", "PValue", "FDR")
]

# Filter Significant Differentially Expressed Genes

sig_EP <- subset(results_EP, FDR < 0.05 & abs(logFC) > 1)

up_EP <- subset(sig_EP, logFC > 1)

down_EP <- subset(sig_EP, logFC < -1)
# Export Results
write.csv(results_EP, "results/edgeR/Early_vs_Proliferating_edgeR_All_Genes.csv", row.names = FALSE)

write.csv(sig_EP, "results/edgeR/Early_vs_Proliferating_edgeR_DEGs.csv", row.names = FALSE)

write.csv(up_EP, "results/edgeR/Early_vs_Proliferating_edgeR_Upregulated.csv", row.names = FALSE)

write.csv(
  down_EP,
  "results/edgeR/Early_vs_Proliferating_edgeR_Downregulated.csv",
  row.names = FALSE
)

# Differential Expression
# Comparison: Deep vs Proliferating

contrast_DP <- makeContrasts(Deep - Proliferating, levels = design)

qlf_DP <- glmQLFTest(fit, contrast = contrast_DP)

results_DP <- topTags(qlf_DP, n = Inf)$table

results_DP$Ensembl_ID <- rownames(results_DP)

gene_ids <- sub("\\..*$", "", results_DP$Ensembl_ID)

results_DP$Gene <- mapIds(
  org.Hs.eg.db,
  keys = gene_ids,
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)

results_DP <- results_DP[
  ,
  c("Gene", "Ensembl_ID", "logFC", "logCPM", "F", "PValue", "FDR")
]

sig_DP <- subset(results_DP, FDR < 0.05 & abs(logFC) > 1)

up_DP <- subset(sig_DP, logFC > 1)

down_DP <- subset(sig_DP, logFC < -1)

write.csv(results_DP, "results/edgeR/Deep_vs_Proliferating_edgeR_All_Genes.csv", row.names = FALSE)

write.csv(sig_DP, "results/edgeR/Deep_vs_Proliferating_edgeR_DEGs.csv", row.names = FALSE)

write.csv(up_DP, "results/edgeR/Deep_vs_Proliferating_edgeR_Upregulated.csv", row.names = FALSE)

write.csv(down_DP, "results/edgeR/Deep_vs_Proliferating_edgeR_Downregulated.csv", row.names = FALSE)

# Differential Expression
# Comparison: Deep vs Early

contrast_DE <- makeContrasts(Deep - Early, levels = design)

qlf_DE <- glmQLFTest(fit, contrast = contrast_DE)

results_DE <- topTags(qlf_DE, n = Inf)$table

results_DE$Ensembl_ID <- rownames(results_DE)

gene_ids <- sub("\\..*$", "", results_DE$Ensembl_ID)

results_DE$Gene <- mapIds(
  org.Hs.eg.db,
  keys = gene_ids,
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)

results_DE <- results_DE[
  ,
  c("Gene", "Ensembl_ID", "logFC", "logCPM", "F", "PValue", "FDR")
]

sig_DE <- subset(results_DE, FDR < 0.05 & abs(logFC) > 1)

up_DE <- subset(sig_DE, logFC > 1)

down_DE <- subset(sig_DE, logFC < -1)

write.csv(results_DE, "results/edgeR/Deep_vs_Early_edgeR_All_Genes.csv", row.names = FALSE)

write.csv(sig_DE, "results/edgeR/Deep_vs_Early_edgeR_DEGs.csv", row.names = FALSE)

write.csv(up_DE, "results/edgeR/Deep_vs_Early_edgeR_Upregulated.csv", row.names = FALSE)

write.csv(down_DE, "results/edgeR/Deep_vs_Early_edgeR_Downregulated.csv", row.names = FALSE)
