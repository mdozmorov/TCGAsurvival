library(dplyr)
library(readr)
clin_xena <- read_csv("data.TCGA/XENA_classification.csv")
clin_p53  <- read_tsv("data.TCGA/BRCA_with_TP53_mutation.tsv")
setdiff(clin_p53$`Case ID`, clin_xena$sampleID)

# Create a vector with indicator 
clin_xena_p53 <- rep("p53wt", nrow(clin_xena))
clin_xena_p53[ clin_xena$sampleID %in% clin_p53$`Case ID` ] <- "p53mut"
# Add P53 status
clin_xena_combined <- clin_xena %>% mutate(TP53_mut_status = clin_xena_p53)

# Pre-load BRCA data, sanity check
setdiff(clin_xena_combined$sampleID, clin$AffyID)

# Make the same order
clin_xena_combined <- clin_xena_combined[ match(clin$AffyID, clin_xena_combined$sampleID), ]
all.equal(clin_xena_combined$sampleID, clin$AffyID)

# Make objects like the original ones
clin_xena <- data.frame(AffyID = clin_xena_combined$sampleID, status = clin_xena_combined$`_OS_IND`, OS = clin_xena_combined$`_OS`)
clin_xena_full <- clin_xena_combined[, !(colnames(clin_xena_combined) %in% c("sampleID"))] %>% as.data.frame()
rownames(clin_xena_full) <- clin_xena_combined$sampleID

# Another sanity check
all.equal(expr$AffyID, rownames(clin_xena_full))

# Save the merged file
write.csv(clin_xena_full, "data.TCGA/BRCA_XENA_clinical.csv")
