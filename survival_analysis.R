source("Supplemental_R_script_1.R")

# Prepare expression data
expr <- mtx$merged.dat[ , 4:ncol(mtx$merged.dat)] %>% as.matrix
# Filter out low expressed genes
# Should be more than 90% of non-zero values
ff <- genefilter::pOverA(p = 0.9, A = 0, na.rm = TRUE) 
expr <- expr[, apply(expr, 2, ff)] 
expr <- data.frame(AffyID = mtx$merged.dat$bcr, expr, stringsAsFactors = FALSE)

# Prepare clinical data
clin <- mtx$merged.dat[, 1:3]
colnames(clin)[1] <- "AffyID"

### Run survival analysis for selected genes
kmplot(expr, clin, event_index=2, time_index=3,  affyid = c("BCL2L11", "MYC", "MYCN"), auto_cutoff="true", transform_to_log2 = TRUE, cancer_type = "BRCA")

### Run survival analysis for all genes
kmplot(expr, clin, event_index=2, time_index=3,  affyid = "", auto_cutoff="true", transform_to_log2 = TRUE)

### Run survival analysis for selected genes in all cancers
selected_genes = "IGFBP3"
# All cancers with RNASeq2 data
cancer_RNASeq2 = c("ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "COADREAD", "DLBC", "ESCA", "GBM", "GBMLGG", "HNSC", "KICH", "KIPAN", "KIRC", "KIRP", "LGG", "LIHC", "LUAD", "LUSC", "MESO", "OV", "PAAD", "PCPG", "PRAD", "READ", "SARC", "SKCM", "STAD", "TGCT", "THCA", "THYM", "UCEC", "UCS"); data.type = "RNASeq2"; type = "" 
#  "LAML" - something wrong with merged data
for (cancer_type in cancer_RNASeq2) {
  print(paste0("Processing cancer ", cancer_type))
  mtx <- load_data(disease = cancer_type, data.type = data.type, type = type, data_dir = data_dir, force_reload = FALSE)
  # Prepare expression data
  expr <- mtx$merged.dat[ , 4:ncol(mtx$merged.dat)] %>% as.matrix
  # Filter out low expressed genes
  # Should be more than 90% of non-zero values
  ff <- genefilter::pOverA(p = 0.9, A = 0, na.rm = TRUE) 
  expr <- expr[, apply(expr, 2, ff)] 
  expr <- data.frame(AffyID = mtx$merged.dat$bcr, expr, stringsAsFactors = FALSE)
  
  # Prepare clinical data
  clin <- mtx$merged.dat[, 1:3]
  colnames(clin)[1] <- "AffyID"
  # Run survival analysis for selected genes
  kmplot(expr, clin, event_index=2, time_index=3,  affyid = selected_genes, auto_cutoff="true", transform_to_log2 = TRUE, cancer_type = cancer_type)
}