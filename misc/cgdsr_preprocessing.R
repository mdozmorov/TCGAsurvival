# http://www.cbioportal.org/cgds_r.jsp
# https://cran.r-project.org/web/packages/cgdsr/vignettes/cgdsr.pdf
# install.packages("cgdsr")
library(cgdsr)
# help('cgdsr')
library(dplyr)
library(annotables)
# Remove non-canonical chromosome names
grch38 <- grch38[ !(grepl("_", grch38$chr) | grepl("GL", grch38$chr)), ]
protein_coding_genes <- grch38[ grch38$entrez != "?" & grch38$biotype == "protein_coding", "symbol"] %>% unlist %>% unique %>% sort

### Standard tutorial
# Create CGDS object
mycgds = CGDS("http://www.cbioportal.org/")
test(mycgds)

# Get available case lists (collection of samples) for a given cancer study
mycancerstudy = getCancerStudies(mycgds)
# Sneak peak into the full list
mycancerstudy[1:5, 1:3] 
View(mycancerstudy)
# Targeted search for studies of interest
mycancerstudy[grepl("neuroblastoma.+TARGET", mycancerstudy$name, ignore.case = TRUE), ]
mycancerstudy[grepl("neuroblastoma.+TARGET", mycancerstudy$name, ignore.case = TRUE), c(1)]
mycancerstudy[grepl("neuroblastoma", mycancerstudy$name, ignore.case = TRUE), c("cancer_study_id", "name")]
#     cancer_study_id                                       name
#        nbl_amc_2012 Neuroblastoma (AMC Amsterdam, Nature 2012)
#      nbl_broad_2013       Neuroblastoma (Broad Institute 2013)
#   nbl_ucologne_2015      Neuroblastoma (Broad, Nat Genet 2013)
# nbl_target_2018_pub     Pediatric Neuroblastoma (TARGET, 2018)

# Select one study
mycancerstudy <- "nbl_target_2018_pub"

# Get all substudies (cases) within a selected study
mycaselist = getCaseLists(mycgds, mycancerstudy)
# Sneak peak into the full list
mycaselist[, c("case_list_id", "case_list_name", "case_list_description")]
View(mycaselist)
#                     case_list_id                                    case_list_name                               case_list_description
#          nbl_target_2018_pub_all                                        All tumors                    All tumor samples (1089 samples)
#    nbl_target_2018_pub_discovery                                  DISCOVERY cohort        All tumors in DISCOVERY cohort (590 samples)
#    nbl_target_2018_pub_sequenced                                  Sequenced Tumors      All (Next-Gen) sequenced samples (372 samples)
#          nbl_target_2018_pub_cna                       Tumor Samples with CNA data               All tumors with CNA data (59 samples)
#         nbl_target_2018_pub_mrna Tumor Samples with mRNA data (Agilent microarray) All samples with mRNA expression data (249 samples)
# nbl_target_2018_pub_rna_seq_mrna            Tumor Samples with mRNA data (RNA Seq) All samples with mRNA expression data (143 samples)
#   nbl_target_2018_pub_validation                                 VALIDATION cohort       All tumors in VALIDATION cohort (499 samples)

# Select one substudy
mycaselist <- "nbl_target_2018_pub_mrna"

# Get available genetic profiles
mygeneticprofile = getGeneticProfiles(mycgds, mycancerstudy)
mygeneticprofile[, c("genetic_profile_id", "genetic_profile_name", "genetic_profile_description")]
View(mygeneticprofile)
#                              genetic_profile_id                         genetic_profile_name                                                                                                                                                                       genetic_profile_description
#                      nbl_target_2018_pub_gistic Putative copy-number alterations from GISTIC Putative copy-number calls on 115 cases determined using GISTIC 2.0. Values: -2 = homozygous deletion; -1 = hemizygous deletion; 0 = neutral / no change; 1 = gain; 2 = high level amplification.
#                nbl_target_2018_pub_rna_seq_mrna               mRNA expression (RNA-Seq RPKM)                                                                                                                                Expression levels for 24946 genes in 144 nbl cases (RNA-Seq RPKM).
# nbl_target_2018_pub_rna_seq_mrna_median_Zscores      mRNA Expression z-Scores (RNA Seq RPKM)                                                                          mRNA z-Scores (RNA Seq RPKM) compared to the expression distribution of each gene tumors that are diploid for this gene.
#        nbl_target_2018_pub_mirna_median_Zscores                 mRNA expression (microarray)                                                                                                                                                    mRNA expression levels (Affymetrix microarray)
#         nbl_target_2018_pub_mrna_median_Zscores        mRNA Expression z-Scores (microarray)                                                                    mRNA z-Scores (Agilent microarray) compared to the expression distribution of each gene tumors that are diploid for this gene.
#                   nbl_target_2018_pub_mutations                            Somatic mutations                                                                                                               Somatic mutation status from whole genome or whole exome sequencing on 313 samples.

# Select one data type
# mygeneticprofile <- "nbl_target_2018_pub_mrna" # No longer works, 04-22-2018 created mtx_nbl_target_2018_pub_mrna.rda
mygeneticprofile <- "nbl_target_2018_pub_mrna_median_Zscores"

# Get data slices for a specified list of genes, genetic profile and case list
selected_genes = c("LONRF2")
myprofiledata <- getProfileData(x = mycgds, genes = c(selected_genes), geneticProfiles = mygeneticprofile, caseList = mycaselist)
View(myprofiledata)

# Get gull data matrix
myprofiledata_full <- list()
# Process all protein-coding genes in chunks of 1000
for (i in 1:(length(protein_coding_genes) %/% 1000)) {
  myprofiledata_chunk <- getProfileData(mycgds, protein_coding_genes[(i * 1000 - 999):min((i * 1000), length(protein_coding_genes))], mygeneticprofile, mycaselist)
  myprofiledata_chunk <- t(myprofiledata_chunk)
  myprofiledata_chunk <- myprofiledata_chunk[complete.cases(myprofiledata_chunk), ]
  myprofiledata_full  <- c(myprofiledata_full, list(myprofiledata_chunk))
}

myprofiledata_full <- do.call("rbind", myprofiledata_full)
myprofiledata_full[1:5, 1:5]
dim(myprofiledata_full)

# Get clinical data for the case list
myclinicaldata = getClinicalData(mycgds, mycaselist)
myclinicaldata[1:5, 1:5]
dim(myclinicaldata)
View(myclinicaldata)
colnames(myclinicaldata)
myclinicaldata <- myclinicaldata[ match(colnames(myprofiledata_full), rownames(myclinicaldata)), ]

all.equal(colnames(myprofiledata_full), rownames(myclinicaldata))

mtx_merged.dat <- data.frame(bcr = colnames(myprofiledata_full), status = ifelse(myclinicaldata$OS_STATUS == "LIVING", 0, 1), OS = myclinicaldata$OS_DAYS, t(myprofiledata_full))
rownames(mtx_merged.dat) <- NULL
mtx_merged.dat[1:5, 1:5]

mtx <- list()
mtx$clinical <- myclinicaldata
mtx$merged.dat <- mtx_merged.dat

# General settings
# Path where the downloaded data is stored
data_dir = "/Users/mdozmorov/Documents/Data/GenomeRunner/TCGAsurvival/data" # Mac
cancer = "nbl_target" 
data.type = "2018_pub"
# type = "mrna" # Doesn't work as of 08-15-2018
type = "mrna_median_Zscores" 
FILE = paste0(data_dir, "/mtx_", cancer, "_", data.type, "_", type, ".rda") # R object with data

save(file = FILE, list = c("mtx")) # Save it
