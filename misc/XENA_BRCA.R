# https://xenabrowser.net/datapages/?cohort=TCGA%20Breast%20Cancer%20(BRCA)
# Pre-load BRCA expression matrix and clinical data, `mtx` object
# Check what we have
expr[1:5, 1:5]
clin[1:5, ]
mtx$clinical[1:5, 1:5]
clin_full <- mtx$clinical

library(jsonlite)
library(readr)

clin_xena <- read_tsv("https://tcga.xenahubs.net/download/TCGA.BRCA.sampleMap/BRCA_clinicalMatrix.gz") # Clinical data
# Keep only tumor samples, and delected columnt
clin_xena <- clin_xena[clin_xena$sample_type == "Primary Tumor", c("sampleID", "PAM50Call_RNAseq", "_OS", "_OS_IND", "_RFS", "_RFS_IND", "_TIME_TO_EVENT", "_TIME_TO_EVENT_UNIT", "breast_carcinoma_estrogen_receptor_status", "breast_carcinoma_progesterone_receptor_status", "breast_carcinoma_surgical_procedure_name", "cytokeratin_immunohistochemistry_staining_method_mcrmtstss_ndctr", "gender", "histological_type", "lab_proc_her2_neu_immunohistochemistry_receptor_status", "margin_status", "menopause_status", "pathologic_M", "pathologic_N", "pathologic_T", "pathologic_stage", "person_neoplasm_cancer_status", "primary_lymph_node_presentation_assessment", "radiation_therapy", "sample_type")]

# # Getting expression data
# mtx_xena <- read_tsv("https://tcga.xenahubs.net/download/TCGA.BRCA.sampleMap/HiSeqV2.gz")
# mtx_xena <- mtx_xena[!(grepl("?", mtx_xena$sample, fixed = TRUE)), ] # Exclude genes without gene name
# mtx_xena <- mtx_xena[ order(mtx_xena$sample), ] # Order by gene names
# genes_xena <- make.names(mtx_xena$sample) # Save gene names
# mtx_xena$sample <- NULL # Remove gene column
# # Matching with clinical annotations
# mtx_xena <- mtx_xena[ , colnames(mtx_xena) %in% clin_xena$sampleID] # Keep only samples with clinical annotation
# mtx_xena <- data.frame(AffyID = colnames(mtx_xena), (2^t(mtx_xena) - 1) ) # Make look like expr
# colnames(mtx_xena)[-1] <- genes_xena # Add column names = gene names
# mtx_xena$AffyID <- substr(mtx_xena$AffyID, 1, 12) # Cut sample names
# common_names <- intersect(mtx_xena$AffyID, expr$AffyID)
# common_genes <- intersect(colnames(mtx_xena), colnames(expr))
# mtx_xena <- mtx_xena[ mtx_xena$AffyID %in% expr$AffyID, ]
# mtx_xena <- mtx_xena[match(expr$AffyID, mtx_xena$AffyID), ]
# mtx_xena <- mtx_xena[, colnames(mtx_xena) %in% colnames(expr)]
# mtx_xena <- mtx_xena[, match(colnames(expr), colnames(mtx_xena))]
# rownames(mtx_xena) <- NULL
# mtx_xena[1:5, 1:5]
# all.equal(mtx_xena, expr, tolerance = 0.002)
# clin_xena$sampleID <- substr(clin_xena$sampleID, 1, 12) # Cut sample names
# clin_xena <- clin_xena[ clin_xena$sampleID %in% mtx_xena$AffyID, ]
# clin_xena <- clin_xena[ match(mtx_xena$AffyID, clin_xena$sampleID), ]
# all.equal(clin_xena$sampleID, clin$AffyID)

# Match with expr
clin_xena$sampleID <- substr(clin_xena$sampleID, 1, 12) # Cut sample names
common_names <- intersect(expr$AffyID, clin_xena$sampleID)
# Attach to expr sampleIDs
clin_xena <- left_join(data.frame(sampleID = expr$AffyID), clin_xena, by = "sampleID")
# Attach race
clin_xena <- left_join(clin_xena, data.frame(sampleID = rownames(clin_full), race = clin_full[ , "race"]), by = "sampleID")
all.equal(expr$AffyID, clin_xena$sampleID)

sum(!is.na(clin_xena$PAM50Call_RNAseq))
# 840
table(clin_xena$PAM50Call_RNAseq[ !is.na(clin_xena$PAM50Call_RNAseq) ])
table(clin_xena[, c("PAM50Call_RNAseq", "race")])
# Basal   Her2   LumA   LumB Normal 
# 139     67    419    192     23
write_csv(clin_xena, "data.TCGA/XENA_classification.csv")
