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

# Select genes of interest
selected_genes = c("TMEM212", "IGFBP3") # BRCA
selected_genes = c("BCL2L11", "MYC", "MYCN") # LUAD
selected_genes = c("LAPTM4B", "PIP5K1C") # HNSC


### Run survival analysis for selected genes
kmplot(expr, clin, event_index=2, time_index=3,  affyid = selected_genes, auto_cutoff="true", transform_to_log2 = TRUE, cancer_type = "BRCA")

### Run survival analysis for all genes
kmplot(expr, clin, event_index=2, time_index=3,  affyid = "", auto_cutoff="true", transform_to_log2 = TRUE)

### Run survival analysis for selected genes in all cancers
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

### Run survival analysis for selected genes, in selected cancer, 
### across all combinations of categories in each clinical annotation
# Full clinical information
clin_full <- mtx$clinical
# For each clinical annotation
for (annotation in clinical_annotations) { 
  # Get the number of patients per category in the current annotation
  annotations <- table(clin_full[, annotation], useNA = "no") 
  # How many categories in the current annotation
  num_of_annot_categories <- length(annotations) 
  # How many categories to select at one time
  for (num_of_selected_categories in 1:num_of_annot_categories) {
    # All combinations of categories, m categories selected at one time
    combination_of_categories <- combn(x = names(annotations), m = num_of_selected_categories)
    # For each combination of categories (column)
    for (combination in 1:ncol(combination_of_categories)) {
      # Select patients annotated with this combination of categories
      patients <- rownames(clin_full)[ clin_full[, annotation] %in% combination_of_categories[, combination]]
      # Get their index in the clin and expr matrixes
      index_patients <- which(clin$AffyID %in% patients)
      # If the number of patients annotated with the combination of categories is large enough, proceed
      if (length(index_patients) > 40) {
        print(paste("Processing annotation:", annotation, 
                    ", categories:", paste(combination_of_categories[, combination], collapse = ","),
                    ", number of patients:", length(index_patients)))
        # Get a subset of clinical information for these patients
        clin_selected <- clin[ index_patients, ]
        # Get a subset of expression information for these patients
        expr_selected <- expr[ index_patients, ]
        # For this subset of expression, filter out low expressed genes
        index_genes <- apply(expr_selected %>% dplyr::select(-AffyID), 2, ff) # index of expression values to keep
        expr_selected <- cbind(expr_selected$AffyID, select(expr_selected, -AffyID)[, index_genes]) # patient IDs and expression values to keep
        # Perform actual survival analysis
        kmplot(expr_selected, clin_selected, event_index=2, time_index=3,  affyid = selected_genes, auto_cutoff="true", transform_to_log2 = TRUE, cancer_type = paste(c(cancer, combination_of_categories[, combination]), collapse = "-"))
      }
    }
  }
}


