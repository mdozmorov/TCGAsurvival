### Pan-cancer correlation analysis, selected genes vs. all others, all cancers (or one), no clinical annotations
# All cancers with RNASeq2 data
cancer_RNASeq2 = c("ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "COADREAD", "DLBC", "ESCA", "GBM", "GBMLGG", "HNSC", "KICH", "KIPAN", "KIRC", "KIRP", "LGG", "LIHC", "LUAD", "LUSC", "MESO", "OV", "PAAD", "PCPG", "PRAD", "READ", "SARC", "SKCM", "STAD", "TGCT", "THCA", "THYM", "UCEC", "UCS"); data.type = "RNASeq2"; type = "" 
# Get correlation matrixes for the gene of interest in each cancer
for (cancer_type in cancer_RNASeq2) {
  print(paste0("Processing cancer ", cancer_type))
  # Prepare expression data
  mtx <- load_data(disease = cancer_type, data.type = data.type, type = type, data_dir = data_dir, force_reload = FALSE)
  expr <- mtx$merged.dat[ , 4:ncol(mtx$merged.dat)] %>% as.matrix
  # Filter out low expressed genes
  # Should be more than 90% of non-zero values
  ff <- genefilter::pOverA(p = 0.9, A = 0, na.rm = TRUE) 
  expr <- expr[, apply(expr, 2, ff)] 
  expr <- data.frame(AffyID = mtx$merged.dat$bcr, expr, stringsAsFactors = FALSE)
  # Correlation analysis
  cor.corr <- apply(expr, 2, function(x) Hmisc::rcorr(as.numeric(x), as.numeric(expr[, which(colnames(expr) == selected_genes)]))[[1]][1, 2])
  cor.pval <- apply(expr, 2, function(x) Hmisc::rcorr(as.numeric(x), as.numeric(expr[, which(colnames(expr) == selected_genes)]))[[3]][1, 2])
  cor.combined <- data.frame(hgnc = names(cor.corr), cor = cor.corr, pval = cor.pval)
  # Save the results
  fileName <- paste0("results/CORR_", selected_genes, ".xlsx")
  if (!file.exists(fileName)) {
    wb <- createWorkbook(fileName)
  } else {
    wb <- loadWorkbook(fileName)
  }
  addWorksheet(wb, sheetName = cancer_type)
  writeData(wb, cor.combined[ order(cor.combined$cor, decreasing = TRUE), ], sheet = cancer_type)
  saveWorkbook(wb, fileName, overwrite = TRUE)
}

# For all cancers only
# Get genes most frequently correlated with the gene of interest (cor > 0.5, pval < 0.05) across all cancers
all_corrs <- list() # List to store cancer-specific correlation matrixes
for (cancer_type in cancer_RNASeq2) { # Go through each file containing cancer-specific correlation matrix
  print(paste0("Processing cancer ", cancer_type))
  mtx <- read.xlsx(fileName, sheet = cancer_type)
  # Get genes positively (cor > 0.5) and significantly (p < 0.05) correlated with the gene of interest
  all_corrs[length(all_corrs) + 1] <- list(mtx[ mtx$cor > 0.5 & mtx$pval < 0.05, grep("^hgnc|^cor|^pval", colnames(mtx))])
}
all_summary <- Reduce(function(...) full_join(..., by = "hgnc"), all_corrs) # Combine all correlation matrixes with sighificant correlations
all_summary <- mutate(all_summary, total=apply(all_summary, 1, function(x) (sum(!is.na(x)) - 1) / 2 )) # Append the total number of cancers having a gene correlated with the gene of interest
all_summary <- all_summary[order(all_summary$total, decreasing = TRUE), grep("^hgnc|^total", colnames(all_summary)) ] # Order by the total
wb <- loadWorkbook(fileName)
addWorksheet(wb, sheetName = "GLOBAL_CORR")
writeData(wb, all_summary, sheet = "GLOBAL_CORR")
saveWorkbook(wb, fileName, overwrite = TRUE)

# Enrichment analysis, all cancers only
library(MDmisc)
library(org.Hs.eg.db)
# Select significant genes
hist(all_summary$total) # Overall distribution
total_cutoff <- as.numeric(quantile(all_summary$total, 0.75)) # Counts of cancers above 75%
sum(all_summary$total >= total_cutoff) # How many genes co-expressed in 75% of all cancers?
# Helper function to save non-empty results
save_res <- function(res, fileName = fileName, wb = wb, sheetName = "KEGG") {
  if (nrow(res) > 0) {
    addWorksheet(wb = wb, sheetName = sheetName)
    writeData(wb, res, sheet = sheetName)
    saveWorkbook(wb, fileName, overwrite = TRUE)
  }
}
# Create (or, load)  Excel file
fileName <- paste0("results/CORR_", selected_genes, ".xlsx")
wb <- loadWorkbook(fileName)
# Gene ontology, molecular function
res <- gene_enrichment(selected = all_summary$hgnc[all_summary$total >= total_cutoff], id="symbol", use="GO", ont="MF")
save_res(res, fileName, wb = wb, sheetName = "GLOBAL_GOMF")
# Gene ontology, biological process 
res <- gene_enrichment(selected = all_summary$hgnc[all_summary$total >= total_cutoff], id="symbol", use="GO", ont="BP")
save_res(res, fileName, wb = wb, sheetName = "GLOBAL_GOBP")
# Gene ontology, cellular component
res <- gene_enrichment(selected = all_summary$hgnc[all_summary$total >= total_cutoff], id="symbol", use="GO", ont="CC")
save_res(res, fileName, wb = wb, sheetName = "GLOBAL_GOCC")
# KEGG canonical pathways
res <- gene_enrichment(selected = all_summary$hgnc[all_summary$total >= total_cutoff], id="symbol", use="KEGG")
save_res(res, fileName, wb = wb, sheetName = "GLOBAL_KEGG")
