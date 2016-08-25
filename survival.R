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
selected_genes = c("TMEM219", "IGFBP3") # BRCA
selected_genes = c("BCL2L11", "MYC", "MYCN") # LUAD
selected_genes = c("LAPTM4B", "PIP5K1C") # HNSC
selected_genes = c("MTDH", "SND1") # LIHC
selected_genes = c("NF1") # LIHC
selected_genes = c("ANLN") # PAN
selected_genes = c("TAF2") # LIHC
selected_genes = c("CPEB2") # BRCA
selected_genes = c("PECAM1", "S1PR1", "SPNS2", "TEK", "TIE1"); cancer_RNASeq2 = c("PAAD")

# View and subset by expression and quantiles of the selected genes
library(ggplot2)
library(reshape2)
selected_genes_expression <- melt(expr[, colnames(expr) %in% selected_genes])
ggplot(selected_genes_expression, aes(x = variable, y = log2(value))) + geom_boxplot()
selected_genes_quantiles_up <- with(selected_genes_expression, tapply(value, variable, quantile, probs = 0.75)) # Select upper quantile!
selected_genes_quantiles_lo <- with(selected_genes_expression, tapply(value, variable, quantile, probs = 0.25)) # Select lower quantile!
# Subset exprs by top expression
selected_index_up <- list() # Collect boolean indexes for each gene
selected_index_lo <- list() # Collect boolean indexes for each gene
for (gene in selected_genes) {
  ind_up <- expr[, gene] > selected_genes_quantiles_up[ gene ] # True, if expressed above the selected upper quantile
  ind_lo <- expr[, gene] < selected_genes_quantiles_lo[ gene ] # True, if expressed below the selected lower quantile
  selected_index_up <- c(selected_index, list(ind_up))
  selected_index_lo <- c(selected_index, list(ind_lo))
}
ind_up <- apply(as.data.frame(selected_index_up), 1, all) # Collapse, all selected genes should be expressed in the upper quantile
ind_lo <- apply(as.data.frame(selected_index_lo), 1, all) # Collapse, all selected genes should be expressed in the lower quantile
ind <- ind_up | ind_lo # One or the other
sum(ind) # How many patients
expr <- expr[ind, ] 
clin <- clin[ind, ]

### Analysis 1: Selected genes, selected cancers, no clinical annotations
kmplot(expr, clin, event_index=2, time_index=3,  affyid = selected_genes, auto_cutoff="true", transform_to_log2 = TRUE, cancer_type = cancer)

### All genes, selected cancers, no clinical annotations
kmplot(expr, clin, event_index=2, time_index=3,  affyid = "", auto_cutoff="true", transform_to_log2 = TRUE)


### Analysis 2: Selected genes, all (or selected) cancers, no clinical annotations
# All cancers with RNASeq2 data
cancer_RNASeq2 = c("ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "COADREAD", "DLBC", "ESCA", "GBM", "GBMLGG", "HNSC", "KICH", "KIPAN", "KIRC", "KIRP", "LGG", "LIHC", "LUAD", "LUSC", "MESO", "OV", "PAAD", "PCPG", "PRAD", "READ", "SARC", "SKCM", "STAD", "TGCT", "THCA", "THYM", "UCEC", "UCS"); data.type = "RNASeq2"; type = "" 
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
### Plot the results of one gene across all cancers
library(dplyr)
library(ggplot2)
# Read in analysis natrix
mtx <- read.table("res/global_stats.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE, fill = TRUE)
# Add -log10-transformed p-value
mtx <- mtx %>% mutate(log10.pval = -1 * log10(p.value))
# Print into PDF
pdf(paste0("res/", selected_genes, "_all_TCGA_cancers.pdf"))
mtx %>% subset(Gene == selected_genes) %>% 
  ggplot(aes(x = factor(Cancer, levels(factor(Cancer))[order(log10.pval)]), y = log10.pval)) + 
  geom_bar(stat = "identity") + 
  theme(legend.position="none") +
  labs(x="Cancer", y="-log10(p-value)") +
  coord_flip()
dev.off()


### Analysis 3: Selected genes, all (or, selected) cancers, all unique categories
for (cancer_type in cancer_RNASeq2) {
  print(paste0("Processing cancer ", cancer_type))
  mtx <- load_data(disease = cancer_type, data.type = data.type, type = type, data_dir = data_dir, force_reload = FALSE)
  clinical_annotations <- summarize_data(mtx = mtx)
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
      # Select patients annotated with this categories
      patients <- rownames(clin_full)[ clin_full[, annotation] %in% names(annotations)[ num_of_selected_categories] ]
      # Get their index in the clin and expr matrixes
      index_patients <- which(clin$AffyID %in% patients)
      # If the number of patients annotated with the combination of categories is large enough, proceed
      if (length(index_patients) > 40) {
        print(paste("Processing annotation:", annotation, 
                    ", categories:", names(annotations)[ num_of_selected_categories],
                    ", number of patients:", length(index_patients)))
        # Get a subset of clinical information for these patients
        clin_selected <- clin[ index_patients, ]
        # Get a subset of expression information for these patients
        expr_selected <- expr[ index_patients, ]
        # For this subset of expression, filter out low expressed genes
        index_genes <- apply(expr_selected %>% dplyr::select(-AffyID), 2, ff) # index of expression values to keep
        expr_selected <- cbind(expr_selected$AffyID, dplyr::select(expr_selected, -AffyID)[, index_genes]) # patient IDs and expression values to keep
        # Perform actual survival analysis
        kmplot(expr_selected, clin_selected, event_index=2, time_index=3,  affyid = selected_genes, auto_cutoff="true", transform_to_log2 = TRUE, cancer_type = paste(c(cancer_type, annotation, names(annotations)[ num_of_selected_categories] ), collapse = "-"))
      }
    }
  }
}


### Selected genes, selected cancers, all combinations of clinical annotations
clin_full <- mtx$clinical # Full clinical information
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
        kmplot(expr_selected, clin_selected, event_index=2, time_index=3,  affyid = selected_genes, auto_cutoff="true", transform_to_log2 = TRUE, cancer_type = paste(c(cancer, annotation, combination_of_categories[, combination]), collapse = "-"))
      }
    }
  }
}


