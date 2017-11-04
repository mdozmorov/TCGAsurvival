## Get clinical information

library(RTCGA)
library(dplyr)

# Which dates (snapshots) are available
checkTCGA('Dates')
# "2016-01-28"
# dir.create( "data2" ) # name of a directory in which data will be stored
releaseDate <- "2016-01-28" # "2015-11-01" # Select one

# Whic cohorts are available?
(cohorts <- infoTCGA() %>% 
    rownames() %>% 
    sub("-counts", "", x=.))
cohorts <- "SARC" # Select one


sapply( cohorts, function(element){
  tryCatch({
    downloadTCGA( cancerTypes = element, 
                  destDir = "data2", 
                  date = releaseDate )},
    error = function(cond){
      cat("Error: Maybe there weren't clinical data for ", element, " cancer.\n")
    }
  )
})

# Shortening paths and directories
list.files( "data2") %>% 
  file.path( "data2", .) %>%
  file.rename( to = substr(.,start=1,stop=50))
# Removing NA files from data2 folder
list.files( "data2") %>%
  file.path( "data2", .) %>%
  sapply(function(x){
    if (x == "data2/NA")
      file.remove(x)      
  })
# Paths to clinical data
cohorts %>%
  sapply(function(z){
    list.files("data2") %>%
      file.path("data2", .) %>%
      grep(paste0("_",z,"\\."), x = ., value = TRUE) %>%
      file.path(., list.files(.)) %>%
      grep("clin.merged.txt", x = ., value = TRUE) %>%
      assign(value = .,
             x = paste0(z, ".clinical.path"),
             envir = .GlobalEnv)
  })
# Reading clinical data using readTCGA
ls() %>%
  grep("clinical\\.path", x = ., value = TRUE) %>% 
  sapply(function(element){
    tryCatch({
      readTCGA(get(element, envir = .GlobalEnv),
               dataType = "clinical") -> clinical_file
      
      ## remove non-ASCII strings:
      for( i in 1:ncol(clinical_file)){
        clinical_file[, i] <- iconv(clinical_file[, i],
                                    "UTF-8", "ASCII", sub="")
      } 
      
      assign(value = clinical_file,
             x = sub("\\.path", "", x = element),
             envir = .GlobalEnv )
    }, error = function(cond){
      cat(element)
    }) 
    invisible(NULL)
  }    
  )

## SARC-specific
SARC.colnames <- data.frame(SARC = colnames(SARC.clinical))
View(SARC.colnames)
table(SARC.clinical$patient.tumor_samples.tumor_sample.tumor_histologies.tumor_histology.histological_type)

colnames(SARC.clinical)[grepl("type", colnames(SARC.clinical))]

## BRCA-specific
BRCA.colnames <- data.frame(BRCA = colnames(BRCA.clinical))
View(BRCA.colnames)

colnames(BRCA.clinical)[grepl("progesterone", colnames(BRCA.clinical)) &
                          grepl("status", colnames(BRCA.clinical))]
colnames(BRCA.clinical)[grepl("receptor_status", colnames(BRCA.clinical)) &
                          grepl("patient.breast_carcinoma", colnames(BRCA.clinical))]

table(BRCA.clinical[, "patient.breast_carcinoma_estrogen_receptor_status"])
table(BRCA.clinical[, "patient.lab_proc_her2_neu_immunohistochemistry_receptor_status"])
table(BRCA.clinical[, "patient.breast_carcinoma_progesterone_receptor_status"])

xtabs(~ patient.breast_carcinoma_estrogen_receptor_status + 
        patient.breast_carcinoma_progesterone_receptor_status +
        patient.lab_proc_her2_neu_immunohistochemistry_receptor_status, data = BRCA.clinical)
