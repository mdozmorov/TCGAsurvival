# Scripts to extract TCGA data for survival analysis.

- `TCGA2SynTarget.R` - prepares TCGA data for SynTarget format. Now, general data processing script

- `Supplemental_R_script_1.R` - a modified script to run gene-specific or global survival analysis, from [http://kmplot.com](http://kmplot.com), [Source](http://kmplot.com/analysis/studies/Supplemental%20R%20script%201.R)

- `survival.R` - an example of running survival analyses

- `correlation.R` - co-expression analysis of selected genes vs. others, in all cancers

- `TCGA_DEGs.Rmd` - differential expression analysis of TCGA cohorts separated into groups with high/low expression of selected genes

# TCGA data

Public data is available through the [TCGA2STAT R package](http://www.liuzlab.org/TCGA2STAT/).

- [Cancer types](http://www.liuzlab.org/TCGA2STAT/CancerDataChecklist.pdf)
- [Data types](http://www.liuzlab.org/TCGA2STAT/DataPlatforms.pdf)
- [Clinical values](http://www.liuzlab.org/TCGA2STAT/ClinicalVariables.pdf)

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


