# Load Packages

library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(enrichplot)
library(ggplot2)

# Create Output Folder

dir.create("results/GO", recursive = TRUE, showWarnings = FALSE)

# GO Function

run_GO <- function(input_file, comparison){

  deg <- read.csv(input_file, stringsAsFactors = FALSE)

  ensembl <- sub("\\..*$", "", deg$Ensembl_ID)

  entrez <- mapIds(
    org.Hs.eg.db,
    keys = ensembl,
    column = "ENTREZID",
    keytype = "ENSEMBL",
    multiVals = "first"
  )

  entrez <- na.omit(entrez)

  go_BP <- enrichGO(
    gene = entrez,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05,
    readable = TRUE
  )

  go_results <- as.data.frame(go_BP)

  write.csv(go_results, file = paste0("results/GO/", comparison, "_GO_BP.csv"), row.names = FALSE)

  pdf(paste0("results/GO/", comparison, "_Dotplot.pdf"), width = 10, height = 8)

  print(dotplot(go_BP, showCategory = 20) + ggtitle(comparison))

  dev.off()

  pdf(paste0("results/GO/", comparison, "_Barplot.pdf"), width = 10, height = 8)

  print(barplot(go_BP, showCategory = 20) + ggtitle(comparison))

  dev.off()

  return(go_BP)

}

# Run GO Analyses

go_EP_up <- run_GO(
  "results/edgeR/Early_vs_Proliferating_edgeR_Upregulated.csv",
  "Early_vs_Proliferating_Up"
)

go_EP_down <- run_GO(
  "results/edgeR/Early_vs_Proliferating_edgeR_Downregulated.csv",
  "Early_vs_Proliferating_Down"
)

go_DP_up <- run_GO(
  "results/edgeR/Deep_vs_Proliferating_edgeR_Upregulated.csv",
  "Deep_vs_Proliferating_Up"
)

go_DP_down <- run_GO(
  "results/edgeR/Deep_vs_Proliferating_edgeR_Downregulated.csv",
  "Deep_vs_Proliferating_Down"
)

go_DE_up <- run_GO("results/edgeR/Deep_vs_Early_edgeR_Upregulated.csv", "Deep_vs_Early_Up")

go_DE_down <- run_GO("results/edgeR/Deep_vs_Early_edgeR_Downregulated.csv", "Deep_vs_Early_Down")
