# Scripts to extract TCGA data for survival analysis.

* [Data Description](#data-description)
* [Data preparation](#data-preparation)
* [Analysis examples](#analysis-examples)
* [Analysis scripts](#analysis-scripts)
* [Misc scripts](#misc-scripts)
* [TCGA data](#tcga-data)
  * [OvarianCancerSubtypes](#OvarianCancerSubtypes)
  * [ProteinAtlas](#ProteinAtlas)
  * [brca_mbcproject_wagle_2017](#brca_mbcproject_wagle_2017)
  * [TCGA_Ovarian](#TCGA_Ovarian)

For more cancer-related notes, see https://github.com/mdozmorov/Cancer_notes

## Data description

Scripts are being transitioned to use [curatedTCGAData](https://bioconductor.org/packages/curatedTCGAData/) and [TCGAutils](https://bioconductor.org/packages/TCGAutils/) R packages.

- Public data is available through the [TCGA2STAT R package](http://www.liuzlab.org/TCGA2STAT/), [GitHub repo](https://github.com/zhandong/TCGA2STAT). First, install `BiocManager::install("CNTools")`, clone the repository `git clone https://github.com/zhandong/TCGA2STAT`, and install from source `install.packages("TCGA2STAT_1.2.tar.gz", repos = NULL, type = "source")`
- [Cancer types](data/CancerDataChecklist_AppendixA.PDF)
- [Data types](data/DataPlatforms_AppendixB.PDF)
- [Clinical values](data/ClinicalVariables_AppendixC.PDF)
- [Number of samples per cancer](data.TCGA/TCGA_cancer_counts.csv)

## Data preparation

First, get the data locally using `misc/TCGA_preprocessing.R` script.

- Create a folder on a local computer
- Change the `data_dir` variable with the path where the downloaded data is stored
- Run the file line-by-line, or source it
- By default, RNA-seq data for all cancers will be downloaded and saved as `*.rda` files
- **In all other scripts, change the `data_dir` variable to the path where the downloaded data is stored**

## Analysis examples

- [TNMplot.Rmd](examples/TNMplot.pdf) - differential gene expression analysis in Tumor, Normal and Metastatic Breast Cancer. Reimplementation of online service [tnmplot.com/](https://tnmplot.com/) by Bartha, Áron, and Balázs Győrffy. “[TNMplot.Com: A Web Tool for the Comparison of Gene Expression in Normal, Tumor and Metastatic Tissues](https://doi.org/10.1101/2020.11.10.376228)” 
- [PAM50_EA_AA.Rmd](examples/PAM50_EA_AA.pdf) - Breast cancer, gene expression analysis in 'black or african-american' and 'white' cohorts, in PAM50 subtypes
- [Survival analysis summary, "survival.Rmd", then "TCGA_summary.Rmd"](examples/TCGA_summary.pdf)
- [Differential expression analysis results, "TCGA_DEGs.Rmd"](examples/TCGA_DEGs_MIA.pdf), [Example Exel output](examples/TCGA_DEGs_MIA.xlsx) 
- [Expression analysis summary, "TCGA_expression.Rmd"](examples/TCGA_expression.pdf)
- [Correlation analysis results, "TCGA_correlations.Rmd"](examples/TCGA_correlations_MIA.pdf), [Example Excel output](examples/TCGA_correlations_MIA.xlsx) 
- [CNV analysis of two genes, survival and differential expression, "TCGA_CNV.Rmd"](examples/TCGA_CNV.pdf)

## Analysis scripts

- **In all other scripts, change Path where the downloaded data is stored, `data_dir` variable**

- `survival.Rmd` - a pipeline to run survival analyses for all cancers. Adjust settings `cancer = "BRCA"` and `selected_genes = "IGFBP3"` to the desired cancer and gene IDs. These IDs should be the same in `TCGA_summary.Rmd` that'll summarize the output into [Survival analysis summary](examples/TCGA_summary.pdf). Note if `subcategories_in_all_cancers <- TRUE`, survival analysis is done for all subcategories and all cancers, time consuming.
    - `Analysis 1` - Selected genes, selected cancers, no clinical annotations. Results are in `<selected_genes>.<cancer>.Analysis1` folder.
    - `Exploratory` - All genes, selected cancers, no clinical annotations. Not run by default.
    - `Analysis 2` - Selected genes, all (or selected) cancers, no clinical annotations. Results are in `<selected_genes>.<cancer>.Analysis2` folder.
    - `Analysis 3` - Selected genes, all (or, selected) cancers, all unique clinical (sub)groups. Results are in `<selected_genes>.<cancer>.Analysis3` folder. Open file `global_stats.txt` in Excel, sort by p-value (log-rank test) and explore in which clinical (sub)groups expression of the selected gene affects survival the most.
    - `Analysis 4` - Selected genes, selected cancers, all combinations of clinical annotations. Not run by default.
    - `Analysis 5` - Analysis 5: Clinical-centric analysis. Selected cancer, selected clinical subcategory, survival difference between all pairs of subcategories. Only run for BRCA and OV cancers. Results are in `<selected_genes>.<cancer>.Analysis5`
    - `Analysis 6` - Dimensionality reduction of a gene signature across all cancers using NMF, PCA, or FA For each cancer, extracts gene expression of a signature, reduces its dimensionality, plots a heatmap sorted by the first component, biplots, saves eigenvectors in files named after cancer, signature, method. They are used in `correlations.Rmd`. Not run by default

- `survival_Neuroblastoma.Rmd` - survival analysis for Neuroblastoma samples from TARGET database. Prepare the data with `misc/cgdsr_preprocessing.R`, see Methods section for data description.

- `TCGA_summary.Rmd` - summary report for the `survival.Rmd` output. In which cancers, and clinical subgroups, expression of the selected gene affects survival the most. Change `cancer = "BRCA"` and `selected_genes = "IGFBP3"` to the desired cancer and gene IDs. Uses results from `<selected_genes>.<cancer>.Analysis*` folders. [Survival analysis summary](examples/TCGA_summary.pdf)

- `TCGA_CNV.Rmd` - Separate samples based on copy number variation of one or several genes, do survival and differential expression analysis on the two groups, and KEGG enrichment. An ad hoc analysis, requires manual intervention.

- `TCGA_stemness.Rmd` - correlation of a selected gene with stemness indices, for details, see Malta, Tathiane M., Artem Sokolov, Andrew J. Gentles, Tomasz Burzykowski, Laila Poisson, John N. Weinstein, Bożena Kamińska, et al. “Machine Learning Identifies Stemness Features Associated with Oncogenic Dedifferentiation.” Cell 173, no. 2 (April 2018): 338-354.e15. https://doi.org/10.1016/j.cell.2018.03.034. [Results example PDF](examples/TCGA_stemness.pdf)

- `TCGA_expression.Rmd` - Expression of selected genes across all TCGA cancers. Used for comparing expression of two or more genes. Change `selected_genes <- "XXXX"`, can be multiple. Generates a PDF file with a barplot of log2-expression of selected genes across all cancers, with standard errors. [Example](examples/TCGA_expression.pdf)

- `TCGA_correlations.Rmd` - Co-expression analysis of selected gene vs. all others, in selected cancers. Genes best correlating with the selected gene may share common functions, described in the KEGG canonical pathway analysis section. Gene counts are converted to TPM. Multiple cancers, with the ComBat batch correction for the cohort effect. Change `selected_genes <- "XXXX"` and `cancer <- "YYYY"` variables. The run saves two RData objects, `data/Expression_YYYY.Rda` and `data/Correlation_XXXX_YYYY.Rda`. This speeds up re-runs with the same settings. The full output is saved in `results/Results_XXXX_YYYY.xlsx`. [Example PDF](examples/TCGA_correlations_MIA.pdf), [Example Excel](examples/Results_MIA_BRCA.xlsx)

- `TCGA_correlations_BRCA.Rmd` - Co-expression analysis of selected gene vs. all others, in BRCA stratified by PAM50 annotations. The full output is saved in `results/Results_XXXX_BRCA_PAM50.xlsx`.

- `correlations_one_vs_one.Rmd` - Co-expression analysis of two genes across all cancers. The knitted HTML contains table with correlation coefficients and p-values.

- `TCGA_DEGs.Rmd` - differential expression analysis of TCGA cohorts separated into groups with high/low expression of selected genes. The results are similar to the `correlation` results, most of the differentially expressed genes are also best correlated with the selected genes. This analysis is to explicitly look at the extremes of the selected gene expression and identify KEGG pathways that may be affected. Change `selected_genes = "XXXX"` and `cancer = "YYYY"`. Manually run through line 254 to see which KEGG pathways are enriched. Then, run the code chunk on line 379 to generate a picture of the selected KEGG pathway, [Example](examples/hsa05217.MIA.png), adjust the `![](hsa0YYYY.XXXX.png)` accordingly. Then, recompile the whole document. [Example PDF](examples/TCGA_DEGs_MIA.pdf), [Example Excel](examples/TCGA_DEGs_MIA.xlsx) 

- `TCGA_DEGs_clin_subcategories.Rmd` - differential expression analysis between pairs of clinical subgroups, e.g., within "race" clinical category pairs of subcategories, e.g., "black or african american" vs. "white" subgroups. Output is saved in one Excel file `CANCERTYPE_DEGs_clin_subcategories.xlsx` with pairs of worksheets, one containing DEGs and another containing enrichment results. Data tables have headers describing individual comparisons and results.

- `PPI_Networks.Rmd` - experimenting with extracting and visualizing data from different PPI databases, for a selected gene.

- `Supplemental_R_script_1.R` - a modified script to run gene-specific or global survival analysis, from [http://kmplot.com](http://kmplot.com), [Source](http://kmplot.com/analysis/studies/Supplemental%20R%20script%201.R)

- `TCPA_correlation.Rmd` - experimenting with TCPA data.

# Legacy analyses

- [CELLX analysis summary, "CELLX_analysis.Rmd"](examples/CELLX_analysis.pdf). Compare with published [Table 1. Comparison of IGFBP3 Expression between Cancers and Normal Tissue](https://www.sciencedirect.com/science/article/pii/S1936523316301371?via%3Dihub#t0005), [Paper](https://www.ncbi.nlm.nih.gov/pubmed/27888710)

## Misc scripts

`misc` folder

- `aracne._networks.R` - experimenting with `aracne.networks` R package, https://www.bioconductor.org/packages/release/data/experiment/html/aracne.networks.html

- `calc_feature_length.R` - get length for gene symbols, resolving aliases

- `calcTPM.R` - function to calculate TPMs from gene counts, from https://github.com/AmyOlex/RNASeqBits/tree/master/R

- `clinical_annotation_merge_BRCA.R` - merging `XENA_classification.csv` and `BRCA_with_TP53_mutation.tsv` into `BRCA_XENA_clinical.csv`

- `featureCounts2TPM.Rmd` - convert featureCounts output to gene symbol-annotated TPMs

- `cgdsr.R` - exploring the Cancer Genomic Data Server, http://www.cbioportal.org/study?id=msk_impact_2017, http://www.cbioportal.org/cgds_r.jsp, https://cran.r-project.org/web/packages/cgdsr/vignettes/cgdsr.pdf

- `cgdsr_preprocessing.R` - preprocessing the data to the format used in the scripts. Currently, processes TARGET Neuroblastoma data

- `overlap_significance.R` - simple example of Fisher's exact test

- `PCA.R` - exercises on dimensionality reduction of gene signatures

- `RTCGA.R` - experimenting with `RTCGA` package, https://bioconductor.org/packages/release/bioc/html/RTCGA.html

- `TCGA_preprocessing.R` - utilities for download and formatting of TCGA data. Use `load_data` and `summarize_data` functions to load cancer-specific expression and clinical data.

- `survplot_0.0.7.tar.gz` - package needed for survival plots

- `XENA_BRCA.R` - Exploring data from Xena UCSC genome browser, https://xenabrowser.net/datapages/?cohort=TCGA%20Breast%20Cancer%20(BRCA)


## TCGA data

`data.TCGA` folder. Some data are absent from the repository because of large size - download through links.

- `BRCA_with_TP53_mutation.tsv` - 355 TCGA samples with TP53 mutations, [Source](https://portal.gdc.cancer.gov/exploration?cases_offset=300&cases_size=100&facetTab=mutations&filters=~%28op~%27and~content~%28~%28op~%27in~content~%28field~%27cases.project.project_id~value~%28~%27TCGA-BRCA%29%29%29~%28op~%27in~content~%28field~%27genes.gene_id~value~%28~%27ENSG00000141510%29%29%29%29%29&searchTableTab=cases)

- `CCLE_Cell_lines_annotations_20181226.txt` - CCLE cell line annotations, from https://portals.broadinstitute.org/ccle/data

- `CCR-13-0583tab1.xlsx` - TNBCtype predictions for 163 primary tumors in TCGA considered to be TNBC, classification into six TNBC subtypes. See http://cbc.mc.vanderbilt.edu/tnbc/index.php for details. "UNC" - unclassified. Supplementary table 1 from Mayer, Ingrid A., Vandana G. Abramson, Brian D. Lehmann, and Jennifer A. Pietenpol. “New Strategies for Triple-Negative Breast Cancer--Deciphering the Heterogeneity.” Clinical Cancer Research: An Official Journal of the American Association for Cancer Research 20, no. 4 (February 15, 2014): 782–90. doi:10.1158/1078-0432.CCR-13-0583.

- `Immune_resistant_program.xlsx` - A gene expression program associated with T cell exclusion and immune evasion. Supplementary Table S4 - genes associated with the immune resistance program, described in Methods. Jerby-Arnon, Livnat, Parin Shah, Michael S. Cuoco, Christopher Rodman, Mei-Ju Su, Johannes C. Melms, Rachel Leeson, et al. “A Cancer Cell Program Promotes T Cell Exclusion and Resistance to Checkpoint Blockade.” Cell 175, no. 4 (November 2018): 984-997.e24. https://doi.org/10.1016/j.cell.2018.09.006.

- `gene_signatures_323.xls` - 323 gene signatures from Fan, Cheng, Aleix Prat, Joel S. Parker, Yufeng Liu, Lisa A. Carey, Melissa A. Troester, and Charles M. Perou. “Building Prognostic Models for Breast Cancer Patients Using Clinical Variables and Hundreds of Gene Expression Signatures.” BMC Medical Genomics 4 (January 9, 2011): 3. https://doi.org/10.1186/1755-8794-4-3.

- `PAM50_classification.txt` - sample classification into PAM50 types

- `patientsAll.tsv` - TCGA sample clinical information, including PAM50, from https://tcia.at/home

- `TCGA_489_UE.k4.txt` - Ovarian cancer classification into four subtypes, from `https://github.com/aedin/OvarianCancerSubtypes/data/23257362`

- `TCGA_cancer_counts.csv` - number of samples per cancer. Created by `misc/TCGA_preprocessing.R`

- `TCGA_cancers.xlsx` - TCGA cancer abbreviations, from http://www.liuzlab.org/TCGA2STAT/CancerDataChecklist.pdf

- `TCGA_genes.txt` - genes measured in TCGA RNA-seq experiments

- `TCGA_isoforms.xlsx` - Isoform switching analysis of TCGA data, tumor vs. normal. Consequences, survival prediction. Using IsoformSwitchAnalyzeR R package. [Supplementary Table 1](https://mcr.aacrjournals.org/highwire/filestream/37855/field_highwire_adjunct_files/4/175983_2_unknown_upload_4059088_hqghxc.xlsx) - gene- and isoforms differentially expressed in all cancers. From Vitting-Seerup, Kristoffer, and Albin Sandelin. “The Landscape of Isoform Switches in Human Cancers.” Molecular Cancer Research 15, no. 9 (September 2017): 1206–20. https://doi.org/10.1158/1541-7786.MCR-16-0459.

- `TCGA_purity.xlsx` - Tumor purity estimates for TCGA samples. Tumor purity estimates according to four methods and the consensus method for all TCGA samples with available data. https://www.nature.com/articles/ncomms9971#supplementary-information. Supplementary Data 1 from Aran, Dvir, Marina Sirota, and Atul J. Butte. “Systematic Pan-Cancer Analysis of Tumour Purity.” Nature Communications 6, no. 1 (December 2015). https://doi.org/10.1038/ncomms9971.

- `TCGA_sample_types.xlsx` - Cancer types and subtypes for all TCGA samples. Includes BRCA subtypes, and subtyping of other cancers, where applicable. PMID: 29625050. [Source](https://www.cell.com/cms/10.1016/j.cell.2018.03.035/attachment/ec810e54-7eb5-43f4-9eaf-3bca356bf347/mmc1.xlsx)

- `TCGA_stemness.xlsx` - Supplementary Table 1 - stemness indices for all TCGA samples. Stemness indices built from various data: mRNAsi - gene expression-based, EREG-miRNAsi - epigenomic- and gene expression-baset, mDNAsi, EREG-mDNAsi - same but methylation-based, DMPsi - differentially methylated probes-based, ENHsi - enhancer-based. Each stemness index (si) ranges from low (zero) to high (one) stemness. From Malta, Tathiane M., Artem Sokolov, Andrew J. Gentles, Tomasz Burzykowski, Laila Poisson, John N. Weinstein, Bożena Kamińska, et al. “Machine Learning Identifies Stemness Features Associated with Oncogenic Dedifferentiation.” Cell 173, no. 2 (April 2018): 338-354.e15. https://doi.org/10.1016/j.cell.2018.03.034. 

- `TCGA.bib` - BibTex of TCGA-related references

- `TCPA_proteins.txt` - List of 224 proteins profiled by RPPA technology. The Cancer Proteome Atlas, [http://tcpaportal.org/tcpa/](http://tcpaportal.org/tcpa/). Data download: [http://tcpaportal.org/tcpa/download.html](http://tcpaportal.org/tcpa/download.html). Paper: [http://cancerres.aacrjournals.org/content/77/21/e51](http://cancerres.aacrjournals.org/content/77/21/e51)

- `XENA_classification.csv` - PAM50 and other clinical data, [Source](https://xenabrowser.net/datapages/?dataset=TCGA.BRCA.sampleMap/BRCA_clinicalMatrix&host=https://tcga.xenahubs.net) 


### OvarianCancerSubtypes

Sample annotations by ovarian cancer subtypes. https://github.com/aedin/OvarianCancerSubtypes

### ProteinAtlas

Uhlen, Mathias, Cheng Zhang, Sunjae Lee, Evelina Sjöstedt, Linn Fagerberg, Gholamreza Bidkhori, Rui Benfeitas, et al. “A Pathology Atlas of the Human Cancer Transcriptome.” Science (New York, N.Y.) 357, no. 6352 (August 18, 2017). doi:10.1126/science.aan2507. http://science.sciencemag.org/content/357/6352/eaan2507

Supplementary material http://science.sciencemag.org/content/suppl/2017/08/16/357.6352.eaan2507.DC1  
- `Table S2` - summary of tissue specific expression for each gene, in normal and cancer tissues.  
- `Table S6` - summary of survival prognostic value, with a simple "favorable/unfavorable" label for each gene. Each worksheet corresponds to a different cancer.  
- `Table S8` - per-gene summary, in which cancers it is prognostic of survival.  

### brca_mbcproject_wagle_2017

https://www.mbcproject.org/

The Metastatic Breast Cancer Project is a patient-driven initiative. This study includes genomic data, patient-reported data (pre-pended as PRD), medical record data (MedR), and pathology report data (PATH). All of the titles and descriptive text for the clinical data elements have been finalized in partnership with numerous patients in the project. As these data were generated in a research, not a clinical, laboratory, they are for research purposes only and cannot be used to inform clinical decision-making. All annotations have been de-identified. More information is available at www.mbcproject.org.

Data download: http://www.cbioportal.org/study?id=brca_mbcproject_wagle_2017#summary. Data includes 78 patients, 103 samples, sample-specific clinical annotations, Putative copy-number from GISTIC, MutSig regions

### TCGA_Ovarian

- Gene expression, methylation, miRNA expression matrices, from Zhang, Shihua, Chun-Chi Liu, Wenyuan Li, Hui Shen, Peter W. Laird, and Xianghong Jasmine Zhou. “Discovery of Multi-Dimensional Modules by Integrative Analysis of Cancer Genomic Data.” Nucleic Acids Research 40, no. 19 (October 2012): 9379–91. https://doi.org/10.1093/nar/gks725. - Integrative analysis of gene expression, metnylation, miRNA expression, using NMF, implemented in Matlab. Supplementary material from https://academic.oup.com/nar/article/40/19/9379/2414808#supplementary-data. 




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