# Bulk RNA-seq differential expression and enrichment

R scripts for differential expression and functional enrichment from a gene-level count matrix. The repository covers DESeq2 and edgeR comparisons, Gene Ontology enrichment, GSEA, and an exploratory WGCNA preprocessing example.

The workflow starts from gene-level read counts. Raw-read QC, trimming, alignment and featureCounts execution are outside the repository's scope. Input data and research results are not included.

## Analysis overview

| Script | Analysis | Main outputs |
| --- | --- | --- |
| `scripts/RNAseq_DESeq2.R` | DESeq2: Early vs Proliferating and Deep vs Proliferating | PCA/MA plots, annotated DEG tables and a combined DEG table |
| `scripts/RNAseq_edgeR.R` | Expression filtering, TMM normalization and quasi-likelihood tests for the above contrasts plus Deep vs Early | All-gene, DEG, upregulated and downregulated tables for each contrast |
| `scripts/go.R` | GO biological-process over-representation analysis of upregulated and downregulated edgeR genes | Enrichment tables, dot plots and bar plots |
| `scripts/gsea.R` | GO biological-process GSEA using all-gene edgeR tables ranked by logFC | Enrichment tables, dot plots and enrichment plots |
| `scripts/wgcna.R` | Expression filtering, normalization, sample clustering and soft-threshold screening | Sample-clustering PDF; soft-threshold results remain in memory |

DESeq2 and edgeR select DEGs using adjusted significance below 0.05 and absolute log2 fold change greater than 1. GO and GSEA consume **edgeR outputs**, not DESeq2 outputs. The WGCNA script does not construct modules or identify hub genes.

## Input format

Place a featureCounts table at `data/gene_counts.txt`, with:

- A `Geneid` column containing Ensembl gene identifiers suitable for human annotation.
- Six annotation columns followed by exactly nine sample count columns.
- Non-negative, unnormalized integer counts.

**DESeq2, edgeR and WGCNA assign groups by column position**, in this exact order:

| Count columns | Group |
| --- | --- |
| 1–3 | Proliferating |
| 4–6 | Early |
| 7–9 | Deep |

Check the column names against the actual sample metadata before running. The scripts do not validate biological sample assignments, and an incorrect column order can silently produce incorrect comparisons. A sample-name-based metadata interface is not implemented.

## Dependencies

Required R packages: DESeq2, edgeR, org.Hs.eg.db, AnnotationDbi, clusterProfiler, enrichplot, ggplot2 and WGCNA.

The previously tested R and package versions are listed in [`tested-package-versions.txt`](tested-package-versions.txt). This records the tested setup; it is not an environment lockfile or installation script. A reproducible environment definition is not yet included.

## Run the analyses

Run from the repository root after supplying the count table and installing the dependencies:

```sh
# Differential expression: these two analyses are independent.
Rscript scripts/RNAseq_DESeq2.R
Rscript scripts/RNAseq_edgeR.R

# Run only after edgeR completes successfully.
Rscript scripts/go.R
Rscript scripts/gsea.R

# Optional, partial exploratory workflow.
Rscript scripts/wgcna.R
```

Tables and explicitly saved plots are written under `results/DESeq2/`, `results/edgeR/`, `results/GO/`, `results/GSEA/` and `results/WGCNA/`. The edgeR MDS/BCV plots do not have an explicit output device; noninteractive execution may create `Rplots.pdf` in the working directory.

No public demonstration dataset or example outputs are included, so a fresh clone cannot reproduce the reported analysis without additional input data.

## Interpretation and current limitations

- **Study design:** Models use the group labels encoded in the scripts. Verify whether the actual experiment requires batch, pairing or other covariates before interpreting results.
- **Annotation:** Confirm the organism, genome annotation and Ensembl identifier format. Annotation uses `org.Hs.eg.db`; version suffixes are removed for mapping, and ambiguous mappings use the first match.
- **GO background:** No explicit gene universe is supplied. Evaluate whether all tested, mappable genes should define the background for the intended analysis.
- **GSEA ranking:** Genes are ranked by logFC; the first occurrence of each duplicated gene symbol is retained. Review this policy and tied ranks. No explicit random seed is set.
- **Empty enrichment:** Plotting is not guarded against empty enrichment results and may fail when no terms pass the cutoffs.
- **WGCNA:** Includes expression preprocessing, sample clustering and soft-threshold screening. Module detection and downstream network analyses are outside the scope of this script.
