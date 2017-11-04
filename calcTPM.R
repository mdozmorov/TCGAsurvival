## Amy Olex
## 1/20/17
## Function to calculate the TPM values for RNA-seq data.
  
calcTPM <- function(data, feature_length){
  
  ##Calculate the RPK value
  RPK <- matrix(0, nrow=dim(data)[1], ncol=dim(data)[2])
  
  for(row in 1:dim(data)[1]){
    for(col in 1:dim(data)[2]){
      RPK[row,col] <- data[row,col]/feature_length$Length[row]
    }
  }
  
  ##Calculate the sums of each column and divide by 1000000
  scale_factor <- colSums(RPK)/1000000
  
  ##Now divide all values in each column by the scaling factor
  TPM <- t(t(RPK)/scale_factor)
  colnames(TPM) <- names(data)
  row.names(TPM) <- row.names(data)
  return(as.data.frame(TPM))
}

# Create feature_length data frame
# Before getting feature length, pre-load `expr` TCGA gene expression matrix
library(biomaRt)
mart <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
# Get genomic coordinates for all gene symbols from `expr`
id <- gsub(".", "-", colnames(expr), fixed = TRUE) # Replace dots by dashes
genes<-getBM(attributes=c('hgnc_symbol','chromosome_name', "start_position", "end_position"),
             filters='hgnc_symbol', values=id, mart=mart, uniqueRows=T)
genes <- genes[ genes$chromosome_name %in% c(as.character(1:22), "X", "Y"), ] # Keep only auto/sex chromosomes
unconverted <- setdiff(colnames(expr), genes$hgnc_symbol) # Gene names/aliases that didn't get converted
# Get conversion for aliases, https://www.biostars.org/p/14971/
# load the annotation database
library(org.Hs.eg.db)
library(DBI)
# use sql to get alias table and gene_info table (contains the symbols)
# first open the database connection
dbCon <- org.Hs.eg_dbconn()
# write your SQL query
sqlQuery <- 'SELECT * FROM alias, gene_info WHERE alias._id == gene_info._id;'
# execute the query on the database
aliasSymbol <- dbGetQuery(dbCon, sqlQuery) # All aliases-to-gene symbols mapping
result <- aliasSymbol[aliasSymbol[, "alias_symbol"] %in% unconverted, colnames(aliasSymbol) %in% c("alias_symbol", "symbol")] # Get the relevant columns
result <- result[!duplicated(result), ] # Remove duplicates
# Get genomic coordinates for all symbols corresponding to aliases
aliases <-getBM(attributes=c('hgnc_symbol','chromosome_name', "start_position", "end_position"),
             filters='hgnc_symbol', values=result$symbol, mart=mart, uniqueRows=T)
aliases <- aliases[ aliases$chromosome_name %in% c(as.character(1:22), "X", "Y"), ] # Keep only auto/sex chromosomes
aliases <- inner_join(result, aliases, by = c("symbol" = "hgnc_symbol")) # Attach aliases
aliases <- aliases[, colnames(aliases) != "symbol"] # Exclude gene symbols, keep aliases and genomic coordinates
colnames(aliases) <- c('hgnc_symbol','chromosome_name', "start_position", "end_position")
# Combine genomic coordinates of genes and aliases
genes <- rbind(genes, aliases)
# Construct feature_length object
feature_length <- data.frame(Symbol = genes$hgnc_symbol, Length = genes$end_position - genes$start_position, stringsAsFactors = FALSE)
feature_length <- feature_length[!duplicated(feature_length), ] # Remove duplicated length
# Aggregate by minimum length
feature_length <- aggregate(feature_length$Length, by = list(feature_length$Symbol), FUN = min)
colnames(feature_length) <- c("Symbol", "Length") # Rename columns after aggregate
feature_length$Symbol <- make.names(feature_length$Symbol) # Replace dashes by dots
# Still, something is unconverted
unconverted <- setdiff(colnames(expr), feature_length$Symbol) # genes that don't have length
unconverted <- gsub(".", "-", unconverted, fixed = TRUE)
unconverted <- as.data.frame(unconverted)

# Save the resulting object
save(list = c("feature_length"), file = "data/feature_length.Rda")
load(file = "data/feature_length.Rda")


# library(annotables)
# genes <- data.frame(Symbol = grch38$symbol, Length = grch38$end - grch38$start, Chr = grch38$chr, stringsAsFactors = FALSE)
# genes <- genes[ genes$Chr%in% c(as.character(1:22), "X", "Y"), ] # Keep only auto/sex chromosomes
# feature_length <- genes[, colnames(genes) %in% c("Symbol", "Length")]
# feature_length <- feature_length[!duplicated(feature_length), ] # Remove duplicated length
# # Aggregate by minimum length
# feature_length <- aggregate(feature_length$Length, by = list(feature_length$Symbol), FUN = min)
# colnames(feature_length) <- c("Symbol", "Length") # Rename columns after aggregate
# feature_length$Symbol <- make.names(feature_length$Symbol)
# 
# unconverted <- setdiff(colnames(expr), feature_length$Symbol) # genes that don't have length
# unconverted[grep("IL", unconverted)]
# 
# # Save the resulting object
# save(list = c("feature_length"), file = "data/feature_length.Rda")
# load(file = "data/feature_length.Rda")
