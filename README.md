# Bulk RNA-seq Time-Course Workflow (QC, DGE, Functional Analysis): ZRN01 Example #1 Analysis Workflow (RM)

Analysis and pipeline by: Ha-Na Shim 

Date: 09/01/2024

For any inquiries or suggestions to the below pipeline, my email is: hshim1@uchicago.edu

## Table of Contents
- [Data & Inputs](#data--inputs)
- [Part 1: Quality control and sample assessment](#part-1-quality-control-and-sample-assessment)
  - [CPM distribution filtering](#cpm-distribution-filtering)
  - [Distribution of counts](#distribution-of-counts)
  - [DESeq2 normalization assessment](#deseq2-normalization-assessment)
  - [Sample-to-sample correlation & dendrogram](#sample-to-sample-correlation--dendrogram)
  - [PCA (± loadings) & Scree](#pca--loadings--scree)
- [Downstream DGE & Functional Analysis (stubs)](#downstream-dge--functional-analysis-stubs)
- [Reproducibility](#reproducibility)

## Tools/packages used

<!-- Horizontal badges -->
[![R][R-badge]][R-url] [![RStudio][RStudio-badge]][RStudio-url] [![pheatmap][pheatmap-badge]][pheatmap-url] [![ggplot2][ggplot2-badge]][ggplot2-url] [![tidyverse][tidyverse-badge]][tidyverse-url] [![nf-core][nfcore-badge]][nfcore-url] [![Snakemake][Snakemake-badge]][Snakemake-url] [![Salmon][Salmon-badge]][Salmon-url] [![Kallisto][Kallisto-badge]][Kallisto-url] [![MultiQC][MultiQC-badge]][MultiQC-url]

<p align="center">

<!-- Bioconductor hexes you provided -->
<a href="https://bioconductor.org/packages/DESeq2">
  <img src="https://raw.githubusercontent.com/Bioconductor/BiocStickers/master/DESeq2/DESeq2.png" height="125" alt="DESeq2 hex">
</a>
<a href="https://bioconductor.org/packages/clusterProfiler">
  <img src="https://raw.githubusercontent.com/Bioconductor/BiocStickers/master/clusterProfiler/clusterProfiler.png" height="125" alt="clusterProfiler hex">
</a>

<!-- Tidyverse family (official hex-stickers repo) -->
<a href="https://ggplot2.tidyverse.org/">
  <img src="https://raw.githubusercontent.com/rstudio/hex-stickers/main/SVG/ggplot2.svg" height="125" alt="ggplot2 hex">
</a>
<a href="https://www.tidyverse.org/">
  <img src="https://raw.githubusercontent.com/rstudio/hex-stickers/main/SVG/tidyverse.svg" height="125" alt="tidyverse hex">
</a>
<a href="https://dplyr.tidyverse.org/">
  <img src="https://raw.githubusercontent.com/rstudio/hex-stickers/main/SVG/dplyr.svg" height="125" alt="dplyr hex">
</a>

<!-- Bioconductor hexes -->
</a>
<a href="https://bioconductor.org/packages/edgeR">
  <img src="https://raw.githubusercontent.com/Bioconductor/BiocStickers/master/edgeR/edgeR.png" height="125" alt="edgeR hex">
</a>

</p>

<!-- Badge image refs -->
[R-badge]: https://img.shields.io/badge/R-276DC3?style=for-the-badge&logo=r&logoColor=white
[RStudio-badge]: https://img.shields.io/badge/RStudio-75AADB?style=for-the-badge&logo=rstudio&logoColor=white
[pheatmap-badge]: https://img.shields.io/cran/v/pheatmap?style=for-the-badge&label=pheatmap
[ggplot2-badge]: https://img.shields.io/cran/v/ggplot2?style=for-the-badge&label=ggplot2
[tidyverse-badge]: https://img.shields.io/cran/v/tidyverse?style=for-the-badge&label=tidyverse
[nfcore-badge]: https://img.shields.io/badge/nf--core-pipelines-00A98F?style=for-the-badge
[Snakemake-badge]: https://img.shields.io/badge/Snakemake-workflows-3277a8?style=for-the-badge
[Salmon-badge]: https://img.shields.io/badge/Salmon-quantification-ff8066?style=for-the-badge
[Kallisto-badge]: https://img.shields.io/badge/kallisto-quantification-555555?style=for-the-badge
[MultiQC-badge]: https://img.shields.io/pypi/v/multiqc?style=for-the-badge&label=MultiQC

<!-- Link refs -->
[R-url]: https://www.r-project.org/
[RStudio-url]: https://posit.co/products/open-source/rstudio/
[pheatmap-url]: https://cran.r-project.org/package=pheatmap
[ggplot2-url]: https://ggplot2.tidyverse.org/
[tidyverse-url]: https://www.tidyverse.org/
[nfcore-url]: https://nf-co.re/
[Snakemake-url]: https://snakemake.readthedocs.io/
[Salmon-url]: https://combine-lab.github.io/salmon/
[Kallisto-url]: https://pachterlab.github.io/kallisto/
[MultiQC-url]: https://multiqc.info/


# Part 1: Quality control and sample clustering

This section of the workflow is primarily focused on assessing the quality of samples and to identify potential outlier samples. Although multiQC or other RNAseq QC pipelines will likely identify problematic samples upstream of generation of raw counts, sample clustering can also provide biologically insightful results.

 [1a. CPM distribution filtering](#1a-cpm-distribution-filtering)
- [1b. Distribution of p-values across DESeq2 contrasts of interest](#1b-distribution-of-p-values-across-deseq2-contrasts-of-interest)
- [1c. log2FoldChange scatterplots to compare shrinkage vs. no shrinkage](#1c-log2foldchange-scatterplots-to-compare-effect-of-deseq2s-apeglm-shrinkage-vs-no-shrinkage)
- [1c. DESeq2 normalization assessment](#1c-deseq2-normalization-assessment)
- [1d. Sample-to-sample correlation matrix heatmap](#1d-sample-to-sample-correlation-matrix-heatmap)
- [1e. PCA plot (± loadings + Scree)](#1e-pca-plot-with-and-withlout-loadings--scree-plot)
- [1f. Lineplot of gene subtypes across all samples](#1f-lineplot-of-gene-subtypes-across-all-samples)

## 1a. CPM distribution filtering

<p align="center">
  <img src="1-figures/quality_control/qc_f1_f2_comb.svg" alt="CPM distribution F1–F2" width="48%"/>
  <img src="1-figures/quality_control/qc_f3_f4_comb.svg" alt="CPM distribution F3–F4" width="48%"/>
</p>

<h3>Parameters used for CPM filtering</h3>
<p>
  <strong>Criteria:</strong> &gt;3 CPM in at least 3 replicates
</p>

<pre><code class="language-r">
total_genes <- nrow(cpm_matrix)
expressed_genes <- rowSums(cpm_matrix > 3) >= 3
n_expressed <- sum(expressed_genes)
</code></pre>

## 1b. Library size comparison across all samples

![Library size boxplot](1-figures/quality_control/lib_size_boxplot.svg)

All samples appear to share similar library sizes, so there are 

## 1b. Distribution of p-values across DESeq2 contrasts of interest

![p-value distributions](1-figures/quality_control/pvalue_distributions.svg)

Interestingly, in the case of this experiment, our condition of interest (genotype) appears to induce a large transcriptional pertubation. ~9000 genes pass our pvalue filter of 0.05. This result is not necessarily problematic since we will apply a log2FoldChange cutoff, which should reduce potential FPs. 

## 1c. log2FoldChange scatterplots to compare effect of DESeq2's apeglm shrinkage vs. no shrinkage


DESeq2 offers a wi

## 1c. DESeq2 normalization assessment

<p align="center">
  <img src="1-figures/quality_control/dispersion_plot.svg" alt="DESeq2 dispersion plot" width="45%"/>
  <img src="1-figures/quality_control/rle_combined_fig.svg" alt="RLE boxplot visualization" width="45%"/>
</p>

## 1d. sample-to-sample correlation matrix heatmap

<p align="center">
  <img src="1-figures/quality_control/sample_correlation_heatmap.svg" alt="Sample-to-sample correlation matrix heatmap" width="45%"/>
  <img src="1-figures/quality_control/sample_dendrogram.svg" alt="Sample dendrogram" width="45%"/>
</p>

## 1e. PCA plot (with and withlout loadings + SCREE plot

<p align="center">
  <img src="1-figures/quality_control/pca_plot_noloadings.svg" alt="PCA plot without loadings" width="30%"/>
  <img src="1-figures/quality_control/pca_plot_loadings.svg" alt="PCA plot with loadings" width="30%"/>
  <img src="1-figures/quality_control/scree_plot.svg" alt="Scree plot" width="30%"/>
</p>

## 1f. Gene subtype lineplot across

![Gene biotype lineplot](1-figures/quality_control/gene_biotype_lineplot.svg)

