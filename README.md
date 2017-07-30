# Scripts to extract TCGA data for survival analysis.

- `CELLX_analysis.Rmd` - Tumor-normal expression of selected gene in all TCGA cancers. In which cancers expression of the selected gene up- or downregulated the most in tumor vs. normal comparison. See the instructions in the document. Change `gene <- "XXXX"` as needed. Output is the HTML file.

- `cgdsr.R` - exploring the Cancer Genomic Data Server, http://www.cbioportal.org/study?id=msk_impact_2017, http://www.cbioportal.org/cgds_r.jsp, https://cran.r-project.org/web/packages/cgdsr/vignettes/cgdsr.pdf

- `survival.Rmd` - a pipeline to run survival analyses. Based on `survival.R`. Change
    - `Analysis 1` - Selected genes, selected cancers, no clinical annotations. Results are in `res.genes.Analysis1` folder.
    - `Exploratory` - All genes, selected cancers, no clinical annotations. Not run by default.
    - `Analysis 2` - Selected genes, all (or selected) cancers, no clinical annotations. Results are in `res.genes.Analysis2` folder.
    - `Analysis 3` - Selected genes, all (or, selected) cancers, all unique clinical (sub)groups. Results are in `res.genes.Analysis3` folder. Open file `global_stats.txt` in Excel, sort by p-value (Cox proportional hazard analysis) and explore in which clinical (sub)groups expression of the selected gene affects survival the most.
    - `Analysis 4` - Selected genes, selected cancers, all combinations of clinical annotations. Not run by default.
    - `Analysis 5` - Clinical-centric analysis. Selected cancer, selected clinical category and two subcategories, survival difference between the two subcategories. Also, saves a plot with boxplots of log2 expression in each clinical group. Not run by default.
    - `Analysis 6` - Dimensionality reduction of a gene signature across all cancers using NMF, PCA, or FA For each cancer, extracts gene expression of a signature, reduces its dimensionality, plots a heatmap sorted by the first component, biplots, saves eigenvectors in files named after cancer, signature, method. They are used in `correlations.Rmd`

- `TCGA_summary.Rmd` - in which cancers, and clinical subgroups, expression of the selected gene affects survival the most. Search and replace the name of the selected gene, and cancer type. Uses results from `res.genes.Analysis2` and `res.genes.Analysis3` folders. Change `gene <- "XXXX"` as needed. Adjust two `![](res.genes.AnalysisX/XXXX.png)` placeholders.

- `coexpression.Rmd` - Expression of selected genes across all TCGA cancers. Change `selected_genes <- "XXXX"`, can be multiple. Generates an HTML file with a barplot of log2-expression of selected genes across all cancers, with standard errors.

- `correlations.Rmd` - Co-expression analysis of selected gene vs. all others, in selected cancers. Genes best correlating with the selected gene may share common functions, described in the KEGG canonical pathway analysis section. Change `selected_genes <- "XXXX"` and `cancer_RNASeq2 <- "YYYY"` variables. The run saves two RData objects, `res/YYYY_expression_RNASeq2_.Rda` and `results/YYYY_correlation_XXXX_RNASeq2_.Rda`. This speeds up re-runs with the same settings. The full output is saved in `results/YYYY_results_XXXX_RNASeq2_.xlsx`

- `correlations_one_vs_one.Rmd` - Co-expression analysis of two genes across all cancers. The knitted HTML contains table with correlation coefficients and p-values.

- `TCGA_DEGs.Rmd` - differential expression analysis of TCGA cohorts separated into groups with high/low expression of selected genes. The results are similar to the `correlation` results, most of the differentially expressed genes are also best correlated with the selected genes. This analysis is to explicitly look at the extremes of the selected gene expression and identify KEGG pathways that may be affected. Change `selected_genes = "XXXX"` and `cancer = "YYYY"`. Manually run through line 254 to see which KEGG pathways are enriched. Then, run the code chunk on line 379 to generate a picture of the selected KEGG pathway, adjust the `![](hsa0YYYY.XXXX.png)` accordingly. Then, recompile the whole document.

- `overlap_significance.R` - simple example of Fisher's exact test

- `PCA.R` - exercises on dimensionality reduction of gene signatures

- `PPI_Networks.Rmd` - experimenting with extracting and visualizing data from different PPI databases, for a selected gene.

- `TCPA_correlation.Rmd` - experimenting with TCPA data.
- `TCPA_proteins.txt` - list of proteins profiled at http://tcpaportal.org/tcpa/. 

- `TCGA_preprocessing.R` - utilities for download and formatting of TCGA data. Use `load_data` and `summarize_data` functions to load cancer-specific expression and clinical data.

- `Supplemental_R_script_1.R` - a modified script to run gene-specific or global survival analysis, from [http://kmplot.com](http://kmplot.com), [Source](http://kmplot.com/analysis/studies/Supplemental%20R%20script%201.R)

# TCGA data

Public data is available through the [TCGA2STAT R package](http://www.liuzlab.org/TCGA2STAT/).

- [Cancer types](http://www.liuzlab.org/TCGA2STAT/CancerDataChecklist.pdf), also saved in [TCGA_cancers.xlsx](data.TCGA/TCGA_cancers.xlsx)
- [Data types](http://www.liuzlab.org/TCGA2STAT/DataPlatforms.pdf)
- [Clinical values](http://www.liuzlab.org/TCGA2STAT/ClinicalVariables.pdf)
- [TCGA_genes.txt](data.TCGA/TCGA_genes.txt) - genes measured in TCGA RNA-seq experiments
- [PAM50 classification](data.TCGA/PAM50_classification.txt) - sample classification into PAM50 types
- [CCR-13-0583tab1.xlsx](data.TCGA/CCR-13-0583tab1.xlsx) - TNBCtype predictions for 163 primary tumors in TCGA considered to be TNBC, classification into six TNBC subtypes. See http://cbc.mc.vanderbilt.edu/tnbc/index.php for details. "UNC" - unclassified. Supplementary table 1 from Mayer, Ingrid A., Vandana G. Abramson, Brian D. Lehmann, and Jennifer A. Pietenpol. “New Strategies for Triple-Negative Breast Cancer--Deciphering the Heterogeneity.” Clinical Cancer Research: An Official Journal of the American Association for Cancer Research 20, no. 4 (February 15, 2014): 782–90. doi:10.1158/1078-0432.CCR-13-0583.
- `EINAV_INTERFERON_SIGNATURE_IN_CANCER.txt` - A gene expression signature found in a subset of cancer patients suggestive of a deregulated immune or inflammatory response. http://software.broadinstitute.org/gsea/msigdb/geneset_page.jsp?geneSetName=EINAV_INTERFERON_SIGNATURE_IN_CANCER

# TCPA data

The Cancer Proteome Atlas, [http://tcpaportal.org/tcpa/](http://tcpaportal.org/tcpa/). Data download: [http://tcpaportal.org/tcpa/download.html](http://tcpaportal.org/tcpa/download.html). 

List of 224 proteins profiled by RPPA  technology: [TCPA_proteins.txt](data.TCGA/TCPA_proteins.txt)

# External tools

- `TNBCtype` tool to classify triple negative breast cancer samples (microarray gene expression) into six subtypes, http://cbc.mc.vanderbilt.edu/tnbc/index.php

- `genefu` R package for PAM50 classification and survival analysis. https://www.bioconductor.org/packages/release/bioc/html/genefu.html



<!--
# [SynTarget](http://www.chemoprofiling.org/cgi-bin/GEO/cancertarget/web_run_CT.V0.S1.pl) tool

SynTarget data format:

1. expression matrix - text file:
 
sample;SampleId1;SampleId2;SampleId3;...
probeID1;37.7;45,5;67.54;...
probeID2;37.7;45,5;67.54;...

2. mapping probes to genes
 
probeID1;geneID1   (preferably official gene symbols or NCBI ENTREZ IDs)
probeID2;geneID2
..

3. Sample Annotation (survival time + clinical variables):
 
!! The fisrt 3 column names are mandatory, should be "sample_id, surv_time, dead_1_alive_0". The others are arbitrary (if available), and should specify clinical variables names (i.e. "stage", "P53_mutation_status")
 
sample_id, surv_time, dead_1_alive_0, stage, P53_mutation_status, cellularity, lymph_nodes_positive, ...
MB-0101,90.4,0, 2, WT, moderate, 0, HT/RT, 2,..
MB-4832,141.7,1, 0, WT, moderate, 0, RT, 1, ..
MB-5119,59.8,0, 0, NA, moderate, 0, NONE, 2,..
MB-0117,8.2,0, 2, WT, moderate, 1, HT/RT, 2, ..
-->