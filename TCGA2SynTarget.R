library(TCGA2STAT)
library(dplyr)

cancer = "BRCA" # Cancer types: http://www.liuzlab.org/TCGA2STAT/CancerDataChecklist.pdf
data.type = "RNASeq" # Data types: http://www.liuzlab.org/TCGA2STAT/DataPlatforms.pdf
type = "RPKM" # Expression types: http://www.liuzlab.org/TCGA2STAT/DataPlatforms.pdf
# Clinical values: http://www.liuzlab.org/TCGA2STAT/ClinicalVariables.pdf

# Get the data
# mtx <- getTCGA(disease = cancer, data.type = data.type, type = type, clinical = TRUE)
# Save to speed up loading next time
# save(file = paste0("data/mtx_", cancer, ".rda", list = c("mtx"))
load(file = paste0("data/mtx_", cancer, ".rda"))

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


## Create sample annotation matrix
mtx.sample <- mtx$merged.dat[, c("bcr", "OS", "status")] # First 3 columns are c("sample_id", "surv_time", "dead_1_alive_0")
colnames(mtx.sample) <- c("sample_id", "surv_time", "dead_1_alive_0")
# Save gzipped matrix
fileName.gz <- gzfile(paste0("results/mtx_", cancer, "_3sample.txt.gz"), "w")
write.table(mtx.sample, fileName.gz, sep = ";", quote = FALSE, row.names = FALSE)
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
