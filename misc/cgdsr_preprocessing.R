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

# Select one study
mycancerstudy <- "nbl_target_2018_pub"

# Get all substudies (cases) within a selected study
mycaselist = getCaseLists(mycgds, mycancerstudy)

# Sneak peak into the full list
mycaselist[1:5, 1:4]
View(mycaselist)
# Select one substudy
mycaselist <- "nbl_target_2018_pub_mrna"

# Get available genetic profiles
mygeneticprofile = getGeneticProfiles(mycgds, mycancerstudy)
View(mygeneticprofile)
# Select one data type
mygeneticprofile <- "nbl_target_2018_pub_mrna"

# Get data slices for a specified list of genes, genetic profile and case list
selected_genes = c("LONRF2") # BRCA
myprofiledata <- getProfileData(mycgds, c(selected_genes), mygeneticprofile, mycaselist)
View(myprofiledata)

myprofiledata_full <- list()
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
cancer = "nbl_target" 
data.type = "2018_pub"
type = "mrna" 
# Path where the downloaded data is stored
data_dir = "/Users/mdozmorov/Documents/Data/GenomeRunner/TCGAsurvival/data" # Mac
FILE = paste0(data_dir, "/mtx_", cancer, "_", data.type, "_", type, ".rda") # R object with data

save(file = FILE, list = c("mtx")) # Save it
