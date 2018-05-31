library(TCGA2STAT)
library(dplyr)
library(knitr)

# Path where the downloaded data is stored
data_dir = "/Users/mdozmorov/Documents/Data/GenomeRunner/TCGAsurvival/data" # Mac
# data_dir = "/Users/mdozmorov/Documents/nobackup" # Mac
# data_dir = "F:/Data/GenomeRunner/TCGAsurvival/data" # Windows
# data_dir = "D:"
# results_dir = "/Users/mdozmorov/Dropbox" # Path where the results are stored
results_dir = data_dir

# Cancer types: http://www.liuzlab.org/TCGA2STAT/CancerDataChecklist.pdf
# Data types: http://www.liuzlab.org/TCGA2STAT/DataPlatforms.pdf
# Expression types: http://www.liuzlab.org/TCGA2STAT/DataPlatforms.pdf
# Clinical values: http://www.liuzlab.org/TCGA2STAT/ClinicalVariables.pdf

# General settings
data.type = "RNASeq2"
type = "" 
data.type = "miRNASeq" # miRNASeq - "count" for raw read counts (default); "rpmmm" for normalized read counts
type = "rpmmm"
data.type = "Mutation"
type = "somatic"
type="all"
data.type = "Methylation"
type = "450K"
clinical = TRUE

# Breast cancer
cancer = "BRCA" 
# Ovarian cancer
cancer = "OV"
# Liver hepatocellular carcinoma
cancer = "LIHC" 
# Head and Neck squamous cell carcinoma
cancer = "HNSC" 
# Sarcoma
cancer = "SARC" 
# Pancreatic cancer
cancer = "PAAD" 
# Lung adenocarcinoma
cancer = "LUAD"
# Lung squamous cell carcinoma
cancer = "LUSC"

# A function to load TCGA data, from remote repository, or a local R object
load_data <- function(disease = cancer, data.type = data.type, type = type, data_dir = data_dir, force_reload = FALSE) {
  FILE = paste0(data_dir, "/mtx_", disease, "_", data.type, "_", type, ".rda") # R object with data
  if (all(file.exists(FILE), !(force_reload))) {
    # If the data has been previously saved, load it
    load(file = FILE)
  } else {
    # If no saved data exists, get it from the remote source
    mtx <- getTCGA(disease = disease, data.type = data.type, type = type, clinical = TRUE)
    save(file = FILE, list = c("mtx")) # Save it
  }
  return(mtx)
}

# A function to get data overview
summarize_data <- function(mtx = mtx) {
  print(paste0("Dimensions of expression matrix, genex X patients: ", paste(dim(mtx$dat), collapse = " ")))
  print(paste0("Dimensions of clinical matrix, patients X parameters: ", paste(dim(mtx$clinical), collapse = " ")))
  print(paste0("Dimensions of merged matrix, patients X parameters + genes: ", paste(dim(mtx$merged.dat), collapse = " ")))
  print("Head of the merged matrix")
  print(mtx$merged.dat[1:5, 1:10])
  print("Head of the clinical matrix")
  print(mtx$clinical[1:5, 1:7])
  print("List of clinical values, and frequency of each variable: ")
  clin_vars <- apply(mtx$clinical, 2, function(x) length(table(x[ !(is.na(x) & x != "" )]))) %>% as.data.frame()
  # Filter clinical variables to have at least 2, but no more than 10 categories,
  # And they are not dates
  clin_vars <- clin_vars[ as.numeric(clin_vars$.) > 1 & as.numeric(clin_vars$.) < 10 & !grepl("years|days|date|vital", rownames(clin_vars), perl = TRUE) , , drop = FALSE]
  print(kable(clin_vars))
  return(rownames(clin_vars))
}

# A function to create expression matrix
make_expression_matrix <- function(mtx = mtx, disease = cancer, data.type = data.type, type = type, results_dir = results_dir) {
  # transposed expression matrix (genes start at column 4) 
  mtx.expression <- mtx$merged.dat[, 4:ncol(mtx$merged.dat) ] %>% t 
  # Set column names as patient IDs
  colnames(mtx.expression) <- mtx$merged.dat$bcr 
  # Set row names as probe IDs
  rownames(mtx.expression) <- colnames(mtx$merged.dat)[ 4:ncol(mtx$merged.dat) ] 
  # Save gzipped matrix
  fileName.gz <- gzfile(paste0(results_dir, "/mtx_", disease, "_", data.type, "_", type, "_1expression.txt.gz"), "w")
  write.table(mtx.expression, fileName.gz, sep = ";", quote = FALSE)
  close(fileName.gz)
}

# A function to create probe ID - gene symbol mapping
make_mapping_matrix <- function(mtx = mtx, disease = cancer, data.type = data.type, type = type, results_dir = results_dir) {
  mtx.mapping <- data.frame(probeID = colnames(mtx$merged.dat)[ 4:ncol(mtx$merged.dat) ],
                            geneID = colnames(mtx$merged.dat)[ 4:ncol(mtx$merged.dat) ])
  # Save gzipped matrix
  fileName.gz <- gzfile(paste0(results_dir, "/mtx_", disease, "_", data.type, "_", type, "_2mapping.txt.gz"), "w")
  write.table(mtx.mapping, fileName.gz, sep = ";", quote = FALSE, col.names = FALSE, row.names = FALSE)
  close(fileName.gz)
}

# A function to create sample annotation matrix
make_annotation_matrix <- function(mtx = mtx, disease = cancer, data.type = data.type, type = type, clinical_annotations = clinical_annotations, results_dir = results_dir) {
  mtx.sample <- mtx$merged.dat[, c("bcr", "OS", "status")] # First 3 columns are c("sample_id", "surv_time", "dead_1_alive_0")
  colnames(mtx.sample) <- c("sample_id", "surv_time", "dead_1_alive_0")
  # Append selected clinical annotations
  mtx.sample <- left_join(mtx.sample, data.frame(sample_id = rownames(mtx$clinical), mtx$clinical[, colnames(mtx$clinical) %in% clinical_annotations], stringsAsFactors = FALSE), by = "sample_id")
  mtx.sample[ is.na(mtx.sample) ] <- "N/A"
  mtx.sample <- mtx.sample[, apply(mtx.sample, 2, function(x) length(x[ x != "N/A"])) > 40] # Keep clinical annitations with at least 40 patients
  # Save gzipped matrix
  fileName.gz <- gzfile(paste0(results_dir, "/mtx_", disease, "_", data.type, "_", type, "_3sample.txt.gz"), "w")
  write.table(mtx.sample, fileName.gz, sep = ";", quote = FALSE, row.names = FALSE)
  close(fileName.gz)
}


# mtx <- load_data(disease = cancer, data.type = data.type, type = type, data_dir = data_dir, force_reload = FALSE)

# clinical_annotations <- summarize_data(mtx = mtx)

# make_expression_matrix(mtx = mtx, disease = cancer, data.type = data.type, type = type, results_dir = results_dir)

# make_mapping_matrix(mtx = mtx, disease = cancer, data.type = data.type, type = type, results_dir = results_dir)

# make_annotation_matrix(mtx = mtx, disease = cancer, data.type = data.type, type = type, clinical_annotations = clinical_annotations, results_dir = results_dir)

# Get data for all cancers
get_data <- function(cancers, data.type = data.type, type = type, data_dir = data_dir, force_reload) {
  for (cancer in cancers) {
    print(paste0("Processing cancer ", cancer))
    mtx <- load_data(disease = cancer, data.type = data.type, type = type, data_dir = data_dir, force_reload)
    clinical_annotations <- summarize_data(mtx = mtx)
    make_expression_matrix(mtx = mtx, disease = cancer, data.type = data.type, type = type, results_dir = results_dir)
    make_mapping_matrix(mtx = mtx, disease = cancer, data.type = data.type, type = type, results_dir = results_dir)
    make_annotation_matrix(mtx = mtx, disease = cancer, data.type = data.type, type = type, clinical_annotations = clinical_annotations, results_dir = results_dir)
    rm(list = c("mtx", "clinical_annotations"))
  }
}

# All cancers with RNASeq2 data
data.type = "RNASeq2"; type = "" 
cancer_TCGA = c("ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "COADREAD", "DLBC", "ESCA", "GBM", "GBMLGG", "HNSC", "KICH", "KIPAN", "KIRC", "KIRP", "LAML", "LGG", "LIHC", "LUAD", "LUSC", "MESO", "OV", "PAAD", "PCPG", "PRAD", "READ", "SARC", "SKCM", "STAD", "TGCT", "THCA", "THYM", "UCEC", "UCS")
# All cancers with miRNAseq data. Uncomment to get miRNAseq data
# data.type = "miRNASeq"; type = "rpmmm"
# cancer_TCGA = c("ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "COADREAD", "DLBC", "ESCA", "GBM", "GBMLGG", "HNSC", "KICH", "KIPAN", "KIRC", "KIRP", "LAML", "LGG", "LIHC", "LUAD", "LUSC", "MESO", "OV", "PAAD", "PCPG", "PRAD", "READ", "SARC", "SKCM", "STAD", "TGCT", "THCA", "THYM", "UCEC", "UCS")

# sink("TCGA2SynTarget.txt", split = FALSE)
get_data(cancers = cancer_TCGA, data.type = data.type, type = type, data_dir = data_dir, force_reload = TRUE)
# sink(type = "message")
# sink()

# Cleanup intermediate files
unlink(paste0(data_dir, "/*.txt.gz"))
