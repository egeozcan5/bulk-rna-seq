# Partial workflow: preprocessing, sample QC and soft-threshold screening only.
# Module detection, trait associations and hub-gene analysis are not implemented.
# Confirm the nine sample labels against the count-column order before running.

# Load Packages
library(WGCNA)
library(edgeR)

# Create Output Directories
output_dir <- "results/WGCNA"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

dir.create(file.path(output_dir, "QC"), showWarnings = FALSE)

# Read Count Matrix

counts <- read.delim(
  "data/gene_counts.txt",
  comment.char = "#",
  stringsAsFactors = FALSE
)

# Prepare Expression Matrix

gene_ids <- counts$Geneid
expr <- counts[, 7:ncol(counts)]
rownames(expr) <- gene_ids

# Filter Lowly Expressed Genes

group <- factor(
  c(
    "Proliferating",
    "Proliferating",
    "Proliferating",
    "Early",
    "Early",
    "Early",
    "Deep",
    "Deep",
    "Deep"
  )
)

dge <- DGEList(counts = expr, group = group)

keep <- filterByExpr(dge)

dge <- dge[keep, , keep.lib.sizes = FALSE]

# Normalize Expression Matrix

dge <- calcNormFactors(dge)

expr_norm <- cpm(dge, log = TRUE, prior.count = 1)

# Check Data Quality

datExpr <- t(expr_norm)

gsg <- goodSamplesGenes(datExpr, verbose = 0)

if (!gsg$allOK) {
  datExpr <- datExpr[
    gsg$goodSamples,
    gsg$goodGenes
  ]
}

# Sample Clustering

sampleTree <- hclust(dist(datExpr), method = "average")

pdf(file.path(output_dir, "QC", "SampleClustering.pdf"), width = 10, height = 8)

plot(
  sampleTree,
  main = "Sample Clustering",
  sub = "",
  xlab = "",
  cex.lab = 1.5,
  cex.axis = 1.5,
  cex.main = 2
)

dev.off()

# Choose Soft Threshold

powers <- 1:20

sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 0)
