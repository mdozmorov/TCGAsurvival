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
# Run survival analysis for selected genes
kmplot(expr, clin, event_index=2, time_index=3,  affyid = c("SND1", "MTDH"), auto_cutoff="true", transform_to_log2 = TRUE)
# Run survival analysis for all genes
kmplot(expr, clin, event_index=2, time_index=3,  affyid = "", auto_cutoff="true", transform_to_log2 = TRUE)

