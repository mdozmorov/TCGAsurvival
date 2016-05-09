library(TCGA2STAT)
library(dplyr)

DIR = "/Users/mikhail/Documents/Data/GenomeRunner/TCGAsurvival/data" # Path where the downloaded data is stored

cancer = "BRCA" # Cancer types: http://www.liuzlab.org/TCGA2STAT/CancerDataChecklist.pdf
cancer = "LIHC" # Liver hepatocellular carcinoma
data.type = "RNASeq" # Data types: http://www.liuzlab.org/TCGA2STAT/DataPlatforms.pdf
type = "RPKM" # Expression types: http://www.liuzlab.org/TCGA2STAT/DataPlatforms.pdf
# Clinical values: http://www.liuzlab.org/TCGA2STAT/ClinicalVariables.pdf

# A function to load TCGA data, from remote repository, or a local R object
load_data <- function(disease = cancer, data.type = data.type, type = type, DIR = DIR, force_reload = FALSE) {
  FILE = paste0(DIR, "/mtx_", disease, ".rda") # R object with data
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

summarize_data <- function(mtx = mtx) {
  # Data exploration
  print(paste0("Dimensions of expression matrix, genex X patients: ", paste(dim(mtx$dat), collapse = " ")))
  print(paste0("Dimensions of clinical matrix, patients X parameters: ", paste(dim(mtx$clinical), collapse = " ")))
  print(paste0("Dimensions of merged matrix, patients X parameters + genes: ", paste(dim(mtx$merged.dat), collapse = " ")))
  print("Head of the merged matrix")
  mtx$merged.dat[1:5, 1:10]
  print("Head of the clinical matrix")
  mtx$clinical[1:5, 1:10]
}


mtx <- load_data(disease = cancer, data.type = data.type, type = type, DIR = DIR, force_reload = FALSE)




# Data exploration
# dim(mtx$dat)
# dim(mtx$clinical)
# dim(mtx$merged.dat)
# mtx$merged.dat[1:5, 1:5]


## Create expression matrix
mtx.expression <- mtx$merged.dat[, 4:ncol(mtx$merged.dat) ] %>% t # transposed expression matrix (genes start at column 4) 
colnames(mtx.expression) <- mtx$merged.dat$bcr # Set column names as patient IDs
rownames(mtx.expression) <- colnames(mtx$merged.dat)[ 4:ncol(mtx$merged.dat) ] # Set row names as probe IDs
# Save gzipped matrix
fileName.gz <- gzfile(paste0("results/mtx_", cancer, "_1expression.txt.gz"), "w")
write.table(mtx.expression, fileName.gz, sep = ";", quote = FALSE)
close(fileName.gz)


## Create probe ID - gene symbol mapping
mtx.mapping <- data.frame(probeID = colnames(mtx$merged.dat)[ 4:ncol(mtx$merged.dat) ],
                          geneID = colnames(mtx$merged.dat)[ 4:ncol(mtx$merged.dat) ])
# Save gzipped matrix
fileName.gz <- gzfile(paste0("results/mtx_", cancer, "_mapping.txt.gz"), "w")
write.table(mtx.mapping, fileName.gz, sep = ";", quote = FALSE, col.names = FALSE, row.names = FALSE)
close(fileName.gz)


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

