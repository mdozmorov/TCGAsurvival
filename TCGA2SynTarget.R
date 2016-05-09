library(TCGA2STAT)
library(dplyr)

data_dir = "/Users/mikhail/Documents/Data/GenomeRunner/TCGAsurvival/data" # Path where the downloaded data is stored
results_dir = "/Users/mikhail/Dropbox" # Path where the results are stored

# Cancer types: http://www.liuzlab.org/TCGA2STAT/CancerDataChecklist.pdf
# Data types: http://www.liuzlab.org/TCGA2STAT/DataPlatforms.pdf
# Expression types: http://www.liuzlab.org/TCGA2STAT/DataPlatforms.pdf
# Clinical values: http://www.liuzlab.org/TCGA2STAT/ClinicalVariables.pdf

# Breast cancer
cancer = "BRCA" 
data.type = "RNASeq"
type = "RPKM" 

# Liver hepatocellular carcinoma
cancer = "LIHC" 
data.type = "RNASeq2"
type = ""

# A function to load TCGA data, from remote repository, or a local R object
load_data <- function(disease = cancer, data.type = data.type, type = type, data_dir = data_dir, force_reload = FALSE) {
  FILE = paste0(data_dir, "/mtx_", disease, "_", data.type, "_", type, ".rda") # R object with data
  if (all(file.exists(FILE), !(force_reload))) {
    # If the data has been previously saved, load it
    load(file = FILE)
  } else {
    # If no saved data exists, get it from the remote source
    mtx <- getTCGA(disease = cancer, data.type = data.type, type = type, clinical = TRUE)
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
  print("List of clinical values: ")
  print(colnames(mtx$clinical))
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




mtx <- load_data(disease = cancer, data.type = data.type, type = type, data_dir = data_dir, force_reload = FALSE)

summarize_data(mtx = mtx)

make_expression_matrix(mtx = mtx, disease = cancer, data.type = data.type, type = type, results_dir = results_dir)

make_mapping_matrix(mtx = mtx, disease = cancer, data.type = data.type, type = type, results_dir = results_dir)





## BRCA-specific investigation of clinical parameters
sink(paste0("results/clinical_", cancer, ".txt"))
clinical_annotations <- c("pathologicstage", "pathologyTstage", "pathologyNstage", "pathologyMstage", "radiationtherapy", "histologicaltype", "race", "ethnicity")
for (annotation in clinical_annotations) {
  print("----------------------------------------------------------------")
  print(annotation)
  print(table(mtx$clinical[, annotation]))
}
sink()


## Create sample annotation matrix
mtx.sample <- mtx$merged.dat[, c("bcr", "OS", "status")] # First 3 columns are c("sample_id", "surv_time", "dead_1_alive_0")
colnames(mtx.sample) <- c("sample_id", "surv_time", "dead_1_alive_0")
# Append selected clinical annotations
mtx.sample <- left_join(mtx.sample, data.frame(bcr = rownames(mtx$clinical), mtx$clinical[, colnames(mtx$clinical) %in% clinical_annotations], stringsAsFactors = FALSE), by = "bcr")
mtx.sample[ is.na(mtx.sample) ] <- "N/A"
# Save gzipped matrix
fileName.gz <- gzfile(paste0("results/mtx_", cancer, "_3sample.txt.gz"), "w")
write.table(mtx.sample, fileName.gz, sep = ";", quote = FALSE, row.names = FALSE)
close(fileName.gz)

