# Library
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)

# Output Folder
output_dir <- "results/GSEA"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# GSEA Function

runGSEA <- function(file, comparison){

  # Read edgeR results
  ep <- read.csv(file)

  # Ranked gene list
  geneList <- ep$logFC
  names(geneList) <- ep$Gene

  geneList <- geneList[!is.na(names(geneList))]
  geneList <- geneList[!duplicated(names(geneList))]
  geneList <- sort(geneList, decreasing = TRUE)

  # Run GSEA
  gseaGO <- gseGO(
    geneList = geneList,
    OrgDb = org.Hs.eg.db,
    keyType = "SYMBOL",
    ont = "BP",
    minGSSize = 10,
    maxGSSize = 500,
    pvalueCutoff = 0.05,
    verbose = FALSE
  )

  # Save table
  write.csv(
    gseaGO@result,
    file.path(output_dir, paste0(comparison, "_GSEA_GO.csv")),
    row.names = FALSE
  )

  # Dotplot
  pdf(file.path(output_dir, paste0(comparison, "_GSEA_Dotplot.pdf")), width = 10, height = 8)

  print(dotplot(gseaGO, showCategory = 20))

  dev.off()

  # Enrichment plot
  pdf(file.path(output_dir, paste0(comparison, "_EnrichmentPlot.pdf")), width = 10, height = 8)

  print(gseaplot2(gseaGO, geneSetID = 1))

  dev.off()

  return(gseaGO)

}

# Run GSEA

gsea_EP <- runGSEA(
  "results/edgeR/Early_vs_Proliferating_edgeR_All_Genes.csv",
  "Early_vs_Proliferating"
)

gsea_DP <- runGSEA(
  "results/edgeR/Deep_vs_Proliferating_edgeR_All_Genes.csv",
  "Deep_vs_Proliferating"
)

gsea_DE <- runGSEA("results/edgeR/Deep_vs_Early_edgeR_All_Genes.csv", "Deep_vs_Early")
