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
mycancerstudy[grepl("brca", mycancerstudy$cancer_study_id, ignore.case = TRUE), c(1)]
View(mycancerstudy[grepl("brca", mycancerstudy$cancer_study_id, ignore.case = TRUE), ])

# "brca_metabric" - 2509 samples, microarray data, with clinical annotations
# "brca_bccrc" - CNV, 65 samples               
# "brca_broad" - SNPs, 103 samples             
# "brca_sanger" - SNPs, 100 samples              
# "brca_tcga_pub2015"         
# "brca_tcga_pub"            
# "brca_tcga"                 
# "brca_bccrc_xenograft_2014" - CNV, 160 samples
# "brca_igr_2015" 

# Select one study
mycancerstudy <- "brca_tcga" #  "brca_tcga"

# Get all substudies (cases) within a selected study
mycaselist = getCaseLists(mycgds, mycancerstudy)
# Sneak peak into the full list
mycaselist[1:5, 1:4]
View(mycaselist)
# Select one substudy
mycaselist <- "brca_tcga_all" # "brca_tcga_rna_seq_v2_mrna"

# Get available genetic profiles
mygeneticprofile = getGeneticProfiles(mycgds, mycancerstudy)
View(mygeneticprofile)
# Select one data type
mygeneticprofile <- "brca_tcga_rna_seq_v2_mrna"

# Get data slices for a specified list of genes, genetic profile and case list
selected_genes = c("NF1") # BRCA
selected_genes = c("PSMA1", "PSMA2", "PSMA3", "PSMA4", "PSMA5", "PSMA6", "PSMA7", "PSMB1", "PSMB2", "PSMB3", "PSMB4", "PSMB5", "PSMB6", "PSMB7", "PSMC1", "PSMC2", "PSMC3", "PSMC4", "PSMC5", "PSMC6", "PSMD1", "PSMD2", "PSMD3", "PSMD4", "PSMD5", "PSMD6", "PSMD7", "PSMD8", "PSMD9", "PSMD10", "PSMD11", "PSMD12", "PSMD13", "PSMD14")
myprofiledata <- getProfileData(mycgds, c(selected_genes), mygeneticprofile, mycaselist)
View(myprofiledata)

col2 <- colorRampPalette(c("blue", "white", "yellow"))
col3 <- colorRampPalette(c("blue", "white", "red"))

h <- pheatmap::pheatmap(t( log2(myprofiledata[complete.cases(myprofiledata), ] + 1) ), color = col3(50), scale = "none", cluster_rows = TRUE, cluster_cols = TRUE, clustering_distance_rows = "correlation", clustering_distance_cols = "correlation", clustering_method = "ward.D2",
                   legend = TRUE, show_colnames = FALSE, treeheight_row = 0, cutree_cols = 8)


# Get clinical data for the case list
myclinicaldata = getClinicalData(mycgds, mycaselist)
View(myclinicaldata)
colnames(myclinicaldata)
# DataExplorer::GenerateReport(myclinicaldata, output_dir = paste0("EDA_", mycaselist, "_", mygeneticprofile))

write.csv(cbind(log2(myprofiledata[ h$tree_col$labels, ] + 1), myclinicaldata[h$tree_col$labels, ]), "brca_tcga_mutations_PSM.csv")

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




