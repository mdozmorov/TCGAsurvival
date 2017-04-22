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
selected_genes = c("MTDH", "SND1", "TAF2") # LIHC
selected_genes = c("NF1") # LIHC
selected_genes = c("ANLN") # PAN
selected_genes = c("TAF2") # LIHC
selected_genes = c("CPEB2") # BRCA
selected_genes = c("PECAM1", "S1PR1", "SPNS2", "TEK", "TIE1"); cancer_RNASeq2 = c("PAAD")
selected_genes = c("hsa.mir.142"); cancer_RNASeq2 = c("BRCA")
selected_genes = "MIA"; cancer_RNASeq2 = c("BRCA")
selected_genes = "DYRK1A"
selected_genes = "DCAF7"
selected_genes = "MBD2" # AML
selected_genes = c("FOXP3") # BRCA
selected_genes = c("TNFRSF1B") # TNFR2, BRCA
selected_genes = c("SDC1") # BRCA

# View and subset by expression and quantiles of the selected genes
library(ggplot2)
library(reshape2)
# Define quantile qutoffs
quantile_up <- 0.51 # 0.75
quantile_lo <- 0.49 # 0.25
selected_genes_expression <- melt(expr[, colnames(expr) %in% selected_genes, drop = FALSE])
ggplot(selected_genes_expression, aes(x = variable, y = log2(value))) + geom_boxplot()
selected_genes_quantiles_up <- with(selected_genes_expression, tapply(value, variable, quantile, probs = quantile_up)) # Select upper quantile!
selected_genes_quantiles_lo <- with(selected_genes_expression, tapply(value, variable, quantile, probs = quantile_lo)) # Select lower quantile!
# Subset exprs by top expression
selected_index_up <- list() # Collect boolean indexes for each gene
selected_index_lo <- list() # Collect boolean indexes for each gene
for (gene in selected_genes) {
  ind_up <- expr[, gene] > selected_genes_quantiles_up[ gene ] # True, if expressed above the selected upper quantile
  ind_lo <- expr[, gene] < selected_genes_quantiles_lo[ gene ] # True, if expressed below the selected lower quantile
  selected_index_up <- c(selected_index_up, list(ind_up))
  selected_index_lo <- c(selected_index_lo, list(ind_lo))
}
ind_up <- apply(as.data.frame(selected_index_up), 1, all) # Collapse indexes from multiple selected genes, all selected genes should be expressed in the upper quantile
ind_lo <- apply(as.data.frame(selected_index_lo), 1, all) # Collapse indexes from multiple selected genes, all selected genes should be expressed in the lower quantile
# For survival analysis, combine samples having upper and lower expression of the selected genes
ind <- ind_up | ind_lo # One or the other
sum(ind) # How many patients
expr <- expr[ind, ] 
clin <- clin[ind, ]
# For differential analysis, create group labels
group <- vector(mode = "numeric", length = nrow(expr)) # Empty bector
group[ind_up] <- 1 # Assign numeric groups
group[ind_lo] <- 2
table(group) # How many patients we have
expr <- expr[group != 0, ] # Remove those that are not in quantiles
group <- group[ group != 0 ]

### Analysis 1: Selected genes, selected cancers, no clinical annotations
system("mkdir res") # Create results folder
kmplot(expr, clin, event_index=2, time_index=3,  affyid = selected_genes, auto_cutoff="true", transform_to_log2 = TRUE, cancer_type = cancer, fileType = "png", use_survminer = FALSE)
# Rename the results folder
system("mv res res.genes.Analysis1")
### Exploratory: All genes, selected cancers, no clinical annotations
# kmplot(expr, clin, event_index=2, time_index=3,  affyid = "", auto_cutoff="true", transform_to_log2 = TRUE)


### Analysis 2: Selected genes, all (or selected) cancers, no clinical annotations
# All cancers with RNASeq2 data
cancer_RNASeq2 = c("ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "COADREAD", "DLBC", "ESCA", "GBM", "HNSC", "KICH", "KIPAN", "KIRC", "KIRP", "LIHC", "LUAD", "LUSC", "MESO", "OV", "PAAD", "PCPG", "PRAD", "READ", "SARC", "SKCM", "STAD", "TGCT", "THCA", "THYM", "UCEC", "UCS") # "GBMLGG", "LGG", 
data.type = "RNASeq2"; type = "" 
system("mkdir res") # Create results folder
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
  kmplot(expr, clin, event_index=2, time_index=3,  affyid = selected_genes, auto_cutoff="true", transform_to_log2 = TRUE, cancer_type = cancer_type, fileType = "png", use_survminer = FALSE)
}
### Plot the results of one gene across all cancers
library(dplyr)
library(ggplot2)
# Read in analysis natrix
mtx <- read.table("res/global_stats.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE, fill = TRUE)
# Add -log10-transformed p-value
mtx <- mtx %>% mutate(log10.pval = -1 * log10(p.value))
# Print into PDF
png(paste0("res/", selected_genes, "_all_TCGA_cancers.png"))
mtx %>% subset(Gene == selected_genes) %>% 
  ggplot(aes(x = factor(Cancer, levels(factor(Cancer))[order(log10.pval)]), y = log10.pval)) + 
  geom_bar(stat = "identity") + 
  theme(legend.position="none") +
  labs(x="Cancer", y="-log10(p-value)") +
  coord_flip()
dev.off()
# Rename the results folder
system("mv res res.genes.Analysis2")


### Analysis 3: Selected genes, all (or, selected) cancers, all unique categories
cancer_RNASeq2 = c("BRCA")
system("mkdir res") # Create results folder
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
        kmplot(expr_selected, clin_selected, event_index=2, time_index=3,  affyid = selected_genes, auto_cutoff="true", transform_to_log2 = TRUE, cancer_type = paste(c(cancer_type, annotation, names(annotations)[ num_of_selected_categories] ), collapse = "-"), fileType = "png", use_survminer = FALSE)
      }
    }
  }
}
# Rename the results folder
system("mv res res.genes.Analysis3")


### Analysis 4: Selected genes, selected cancers, all combinations of clinical annotations
system("mkdir res") # Create results folder
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
        kmplot(expr_selected, clin_selected, event_index=2, time_index=3,  affyid = selected_genes, auto_cutoff="true", transform_to_log2 = TRUE, cancer_type = paste(c(cancer, annotation, combination_of_categories[, combination]), collapse = "-"), fileType = "png", use_survminer = FALSE)
      }
    }
  }
}
# Rename the results folder
system("mv res res.genes.Analysis4")

## Plot boxplot of gene expression split by high/low expression
# Make sure correct mtx is loaded
# Check cancer value
# Check selected_gene value
# Get full clinical data
clin.full <- mtx$clinical
# Select major category
clinical_annotations
clinical_category <- "radiationtherapy"
table(clin.full[, clinical_category])
clinical_subcategory <- "yes"
# Cutoff value
cutoff_value <- 10.65 # log2 expression value affecting the survival. from global_stats. DCAF7: BRCA, radiationtherapy-yes - 10.65; :UAD, radiationtherapy-yes - 10.4, 
# Select subset of expression data
selected_samples <- rownames(mtx$clinical)[mtx$clinical[, clinical_category] %in% clinical_subcategory] # Sample names from clinical annotations
mtx_subset <- expr[ mtx$merged.dat$bcr %in% selected_samples, selected_genes] # Expression subset 
mtx_subset <- data.frame(gene = log2(mtx_subset + 1), group = ifelse(log2(mtx_subset + 1) >= cutoff_value, "High", "Low"), stringsAsFactors = FALSE)
# Boxplots split by expression
res <- ggplot(mtx_subset, aes(x = group, y = gene)) + 
  geom_boxplot(aes(fill = group)) +
  xlab("Group") +
  ylab("log2-transformed expression") +
  theme_bw() +
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 14, face = "bold"))
# scale_fill_manual(name = "", values = c("red", "blue"))
plot(res)
ggsave(filename = paste0("res/boxplot_", selected_genes, "_", cancer, "-", clinical_category, "-", clinical_subcategory, ".pdf"), plot = res, width = 3.5, height = 4.5)



