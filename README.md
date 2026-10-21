# Bulk RNA-seq analysis

R scripts for differential expression and enrichment analysis. The contributor confirmed permission to share these scripts. No data or results are included.

## Included work

- DESeq2: two comparisons against Proliferating, PCA/MA plots, annotated DEGs and a combined DEG table (`scripts/RNAseq_DESeq2.R`).
- edgeR: expression filtering, TMM normalization, MDS/BCV plots and quasi-likelihood tests for three supplied contrasts.
- GO: biological-process enrichment of upregulated and downregulated genes.
- GSEA: biological-process enrichment ranked by log fold change.
- WGCNA: partial preprocessing, sample clustering and soft-threshold screening. No modules or hub genes are computed by the supplied script.

QC, trimming, alignment, featureCounts execution scripts were not supplied and are not represented as completed work here.

## Inputs and execution

Run from this directory. Supply a featureCounts table at `data/gene_counts.txt`, with six annotation columns followed by nine count columns. Both edgeR and WGCNA currently assume three Proliferating, three Early and three Deep samples, in that exact order. Verify this against the real sample sheet before execution; a named sample sheet should replace positional labels before release.

```sh
Rscript scripts/RNAseq_DESeq2.R
Rscript scripts/RNAseq_edgeR.R
Rscript scripts/go.R
Rscript scripts/gsea.R
Rscript scripts/wgcna.R
```

edgeR tables are written to `results/edgeR/`; downstream scripts read that same folder. Data and outputs are excluded from Git.

## Dependencies

The supplied scripts reference DESeq2, edgeR, org.Hs.eg.db, AnnotationDbi, clusterProfiler, enrichplot, ggplot2 and WGCNA. Versions used for the local execution test are recorded in `tested-package-versions.txt`. This version list is not an environment lockfile.

## Review before interpreting results

- Confirm human annotation and the Ensembl identifier format against the reference used for counting.
- GO currently uses the package default background. Review whether the background should instead be all tested, mappable genes; this draft preserves the original analysis choice.
- GSEA ranks by logFC and retains the first occurrence of duplicated gene symbols. Review this policy and tied ranks before interpreting enrichment.
- Confirm whether the nine-sample WGCNA exploratory analysis is suitable for the intended scientific use. The supplied workflow is incomplete.
- Verify sample assignments and retain any required authorship acknowledgements. No reuse license has been selected.

## Changes in this draft

Simplified comments and short function calls; removed console-only inspection lines, unused package imports and unused WGCNA output folders. Aligned edgeR output paths with downstream inputs and created required output folders. Statistical thresholds, contrasts, normalization, gene filtering and ranking/background choices were preserved. Reusable GO/GSEA functions and the original WGCNA quality filter remain because they serve the analysis.

Additional empty-enrichment guards were removed at the user's request. Plotting may fail when no terms pass the cutoffs.

All four scripts were run locally against the supplied count table on 21 October 2026. edgeR produced 12 CSV files, GO produced 18 table/plot files, GSEA produced 9 table/plot files, and the partial WGCNA script completed its soft-threshold scan and produced a sample-clustering PDF. No execution-stopping errors were recorded.

Non-fatal warnings included an R/package build-version mismatch, a GO plotting argument warning, tied GSEA ranks and p-value estimation warnings, and WGCNA running sequentially without a parallel backend. Successful execution does not validate sample labels, statistical assumptions or biological conclusions. Results were not included in this repository.

The consolidated DESeq2 script was also executed successfully against the supplied count table. It produced seven CSV tables and three PDFs under `results/DESeq2/`. The original two contrasts, default DESeq2 settings and significance thresholds were retained; redundant annotation and unused normalized-count extraction were removed. GO/GSEA scripts continue to use edgeR outputs.
