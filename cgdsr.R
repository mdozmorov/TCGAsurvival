# http://www.cbioportal.org/cgds_r.jsp
# https://cran.r-project.org/web/packages/cgdsr/vignettes/cgdsr.pdf
# install.packages("cgdsr")
library(cgdsr)
# help('cgdsr')

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
mycancerstudy[grepl("brca", mycancerstudy$cancer_study_id, ignore.case = TRUE), c(1, 2)]
mycancerstudy[grepl("brca", mycancerstudy$cancer_study_id, ignore.case = TRUE), c(3)]

# Select one study
mycancerstudy <- "brca_tcga_pub2015"

# Get all substudies (cases) within a selected study
mycaselist = getCaseLists(mycgds, mycancerstudy)
# Sneak peak into the full list
mycaselist[1:5, 1:4]
View(mycaselist)
# Select one substudy
mycaselist <- "brca_tcga_pub2015_rna_seq_v2_mrna"

# Get available genetic profiles
mygeneticprofile = getGeneticProfiles(mycgds, mycancerstudy)
View(mygeneticprofile)
# Select one data type
mygeneticprofile <- "brca_tcga_pub2015_rna_seq_v2_mrna"

# Get data slices for a specified list of genes, genetic profile and case list
selected_genes = c("SDC1") # BRCA
myprofiledata <- getProfileData(mycgds, c(selected_genes), mygeneticprofile, mycaselist)
View(myprofiledata)

# Get clinical data for the case list
myclinicaldata = getClinicalData(mycgds, mycaselist)
View(myclinicaldata)
colnames(myclinicaldata)
# DataExplorer::GenerateReport(myclinicaldata, output_dir = paste0("EDA_", mycaselist, "_", mygeneticprofile))

## Plot barplot of the selected gene in different subcategories
# Align expression and clinical data
clin <- data.frame(ID = rownames(myclinicaldata), myclinicaldata)
expr <- data.frame(ID = rownames(myprofiledata), myprofiledata)
expr <- expr[ match(clin$ID, expr$ID), ]
all.equal(clin$ID, expr$ID)
# Select clinical category and subcategories
clinical_annotations_selected <- "AJCC_METASTASIS_PATHOLOGIC_PM" # "pathologyMstage"
print(paste0("Number of patients in each subcategory, in the ", clinical_annotations_selected, " category"))
table(clin[, clinical_annotations_selected]) 
# A matrix to plot
mtx_to_plot <- data.frame(Gene = log2(expr[, selected_genes] + 1), Clinical = clin[, clinical_annotations_selected])
# Plot and save
p <- ggplot(melt(mtx_to_plot, id.vars = "Clinical"), aes(x = Clinical, y = value, fill = Clinical)) +
  geom_boxplot() +
  ylab("log2 expression")
plot(p)
ggsave(filename = paste0("res/", cancer, "_", selected_genes, "_", clinical_annotations_selected, ".pdf"), p, width = 5, height = 5)


# Example 1: Association of NF1 copy number alteration and mRNA expression in glioblastoma

df = getProfileData(mycgds, "NF1", c("gbm_tcga_gistic","gbm_tcga_mrna"), "gbm_tcga_all")
head(df)

boxplot(df[,2] ~ df[,1], main="NF1 : CNA status vs mRNA expression", xlab="CNA status", ylab="mRNA expression", outpch = NA)
stripchart(df[,2] ~ df[,1], vertical=T, add=T, method="jitter",pch=1,col='red')

plot(mycgds, "gbm_tcga", "NF1", c("gbm_tcga_gistic","gbm_tcga_mrna"), "gbm_tcga_all", skin = 'disc_cont')

# Example 2: MDM2 and MDM4 mRNA expression levels in glioblastoma

df = getProfileData(mycgds, c("MDM2","MDM4"), "gbm_tcga_mrna", "gbm_tcga_all")
head(df)

plot(df, main="MDM2 and MDM4 mRNA expression", xlab="MDM2 mRNA expression", ylab="MDM4 mRNA expression")

plot(mycgds, "gbm_tcga", c("MDM2","MDM4"), "gbm_tcga_mrna" ,"gbm_tcga_all")

#  Example 3: Comparing expression of PTEN in primary and metastatic prostate cancer tumors

df.pri = getProfileData(mycgds, "PTEN", "prad_mskcc_mrna", "prad_mskcc_primary")
head(df.pri)

df.met = getProfileData(mycgds, "PTEN", "prad_mskcc_mrna", "prad_mskcc_mets")
head(df.met)

boxplot(list(t(df.pri),t(df.met)), main="PTEN expression in primary and metastatic tumors", xlab="Tumor type", ylab="PTEN mRNA expression",names=c('primary','metastatic'), outpch = NA)
stripchart(list(t(df.pri),t(df.met)), vertical=T, add=T, method="jitter",pch=1,col='red')




