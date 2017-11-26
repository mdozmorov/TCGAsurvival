## Amy Olex
## 1/20/17
## Function to calculate the TPM values for RNA-seq data.
#' @param data a data frame with gene IDs rownames and sample IDs colnames.
#' Each cell contains counts
#' @param feature_length a data.frame with "Geneid" and "Length" columns. "Geneid" should correspond to
#' and be the same length and order as rownames(data) 

  
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

