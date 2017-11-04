## Interferon signature
signature <- readLines("/Users/mdozmorov/Documents/Work/VCU_work/Paula_Bos/Interferon/EINAV_INTERFERON_SIGNATURE_IN_CANCER.txt") # Interferon signature
cancer_RNASeq2 <- cancer <- "BRCA"
selected_genes <- "interferon_signature"
signature1 <- signature[ !(signature %in% c("MORC3", "BRD3", "RXRA")) ]
## PSM signature
signature <- openxlsx::read.xlsx("/Users/mdozmorov/Documents/Work/VCU_work/misc/Senthil/PSM_list.xlsx", colNames = FALSE)
signature <- toupper(unlist(signature))
cancer_RNASeq2 <- cancer <- "BRCA"
selected_genes <- "psm_signature"
## FOXP3 signature
signature <- readLines("/Users/mdozmorov/Documents/Work/VCU_work/Paula_Bos/Interferon/FOXP3_signature.txt") # Interferon signature
cancer_RNASeq2 <- cancer <- "BRCA"
selected_genes <- "FOXP3_signature"
signature2 <- signature[ !(signature %in% c("AREG", "HAVCR1", "LAG3")) ]


# Pre-load a matrix with BRCA gene expression in `survival.Rmd`
mtx_selected <- t(log2(mtx$merged.dat[, signature] + 1))
colnames(mtx_selected) <- mtx$merged.dat$bcr
library(pheatmap)

## NMF
method <- "NMF"
library(NMF)
mtx_nmf <- nmf(t(mtx_selected), 3) # with three components
dim(mtx_nmf@fit@W)
dim(mtx_nmf@fit@H)
# Visualize NMF meta-genes
library(reshape2)
mtx_to_plot <- melt(t(mtx_nmf@fit@H))
mtx_to_plot$Var2 <- factor(mtx_to_plot$Var2)
table(mtx_to_plot$Var2)
ggplot(mtx_to_plot, aes(x = Var1, y = value, group = Var2, color = Var2)) + geom_line()
pheatmap(mtx_selected[, order(mtx_nmf@fit@W[, 1])], cluster_cols = FALSE, treeheight_row = FALSE, treeheight_col = FALSE, scale = "row")
mtx_reduced <- mtx_nmf@fit@W
save(mtx_reduced, file = paste0("data/", cancer_RNASeq2, "_", selected_genes, "_", method, ".Rda"))

pheatmap(rbind(log2(mtx$merged.dat[order(mtx_nmf@fit@W[, 1]), "IFI44L"] + 1),
               log2(mtx$merged.dat[order(mtx_nmf@fit@W[, 1]), "FOXP3"] + 1),
               mtx_nmf@fit@W[order(mtx_nmf@fit@W[, 1]), 1]), cluster_rows = FALSE, cluster_cols = FALSE, treeheight_row = FALSE, treeheight_col = FALSE, scale = "row")
cor(cbind(log2(mtx$merged.dat[order(mtx_nmf@fit@W[, 1]), "IFI44L"] + 1),
          log2(mtx$merged.dat[order(mtx_nmf@fit@W[, 1]), "FOXP3"] + 1),
          mtx_nmf@fit@W[order(mtx_nmf@fit@W[, 1]), 1]))

## PCA
# https://cran.r-project.org/web/packages/ggfortify/vignettes/plot_pca.html
method <- "PCA"
library(ggfortify)
autoplot(prcomp(t(mtx_selected), scale. = TRUE, center = TRUE), loadings = TRUE, loadings.label = TRUE)

# What happens when sorting by the first component
mtx_pca <- prcomp(t(mtx_selected), scale. = TRUE, center = TRUE)
summary(mtx_pca)$importance
mtx_pca$rotation[ order(mtx_pca$rotation[, 1], decreasing = TRUE), 1:2]
biplot(mtx_pca)
pheatmap(mtx_selected[, order(mtx_pca$x[, "PC1"])], cluster_cols = FALSE, treeheight_row = FALSE, treeheight_col = FALSE, scale = "row")
# Plotting the expression profile
mtx_pca_x <- mtx_pca$x[ order(mtx_pca$x[, "PC1"], decreasing = TRUE ), ]
rownames(mtx_pca_x) <- seq(1:nrow(mtx_pca$x))
mtx_to_plot <- melt(mtx_pca_x[, 1:3])
mtx_to_plot$Var2 <- factor(mtx_to_plot$Var2)
table(mtx_to_plot$Var2)
ggplot(mtx_to_plot, aes(x = Var1, y = value, color = Var2)) + geom_line()
mtx_reduced <- mtx_pca_x
save(mtx_reduced, file = paste0("data/", cancer_RNASeq2, "_", selected_genes, "_", method, ".Rda"))

## Factor analysis
mtx_fa <- factanal(t(mtx_selected), factors = 3, scores = "regression")
cor(mtx_fa$scores[, 1:3], mtx_pca$x[, 1:3])
autoplot(mtx_fa,  loadings = TRUE, loadings.label = TRUE)


# Find clinical annotations associated with the first eigenvector
# A function to pull out p-value of LM. https://stackoverflow.com/questions/5587676/pull-out-p-values-and-r-squared-from-a-linear-regression
lmp <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
# Make sure the right mtx is preloaded
clinical_annotations <- summarize_data(mtx = mtx) # names of usable clinical annotations
clin <- mtx$clinical[colnames(mtx_selected), clinical_annotations] # Subset clinical annotations by selected samples and annotations
# Vectors to store association results
assoc_radj <- c()
assoc_pval <- c()
# selected_data <- mtx_nmf@fit@W[, 1] # Use NMF
selected_data <- mtx_pca$x[, 1]     # Use PCA
# Test each annotation for association
for (covariate in clinical_annotations){
  pca.lm <- lm( as.numeric(PC1) ~ as.factor(eval(parse(text = covariate))), data = data.frame(PC1 = selected_data, clin ))
  # print(paste(covariate, "accounts for", signif(summary(pca.lm)$adj.r.squared, 5), "variability explained by PC1, p-value", signif(lmp(pca.lm), 5)))
  assoc_radj <- c(assoc_radj, formatC(summary(pca.lm)$adj.r.squared, format = "g", digits = 3))
  assoc_pval <- c(assoc_pval, formatC(lmp(pca.lm), format = "e", digits = 3))
  # pca.lm <- lm( as.numeric(PC2) ~ factor(eval(parse(text = covariate))), data = cbind(sample_annotation, pca$x))
  # print(paste(covariate, "accounts for", signif(summary(pca.lm)$adj.r.squared, 5), "variability explained by the second principle component, p-value", signif(lmp(pca.lm), 5)))
  # pca.lm <- lm( as.numeric(PC3) ~ factor(eval(parse(text = covariate))), data = cbind(sample_annotation, pca$x))
  # print(paste(covariate, "accounts for", signif(summary(pca.lm)$adj.r.squared, 5), "variability explained by the third principle component, p-value", signif(lmp(pca.lm), 5)))
}
assoc_summary <- data.frame(covariate = clinical_annotations, variability = assoc_radj, pval = assoc_pval)
DT::datatable(assoc_summary[order(assoc_summary$variability, decreasing = TRUE), ])

selected_clinical <- "histologicaltype"
mtx_to_plot <- data.frame(Clinical = clin[, selected_clinical], PC1 = selected_data)
mtx_to_plot <- melt(mtx_to_plot)
  
ggplot(mtx_to_plot, aes(x = Clinical, y = value, fill = Clinical)) +
  geom_boxplot() +
  ylab("PC1 relative level") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# Special analysis for PAM50 signature
pam50 <- read.table("/Users/mdozmorov/Documents/Work/GenomeRunner/TCGAsurvival/data.TCGA/PAM50_classification.txt", sep = "\t", header = TRUE)
pam50$Particpant.ID <- sapply(pam50$Particpant.ID, substr, start = 1, stop = 12)
selected_patients <- intersect(pam50$Particpant.ID, names(selected_data))
pam50 <- pam50[ pam50$Particpant.ID %in% selected_patients, ]
mtx_to_plot <- data.frame(Clinical = pam50[ match(selected_patients, pam50$Particpant.ID), "Subtype"], PC1 = selected_data[selected_patients])


## How two signatures relate?
library(FactoMineR)
# Get expression for the first signature
mtx_selected1 <- t(log2(mtx$merged.dat[, signature1] + 1))
colnames(mtx_selected1) <- mtx$merged.dat$bcr
# Get expression for the second signature
mtx_selected2 <- t(log2(mtx$merged.dat[, signature2] + 1))
colnames(mtx_selected2) <- mtx$merged.dat$bcr
# Join them
mtx_selected <- rbind(mtx_selected1, mtx_selected2)
sd_cutoff <- quantile(apply(mtx_selected, 2, sd), p = 0.95) # Select the most variable genes
mtx_selected <- mtx_selected[, apply(mtx_selected, 2, sd) > sd_cutoff] # here
# Append qualitative column
mtx_selected <- data.frame(mtx_selected, signature = c(rep("IFN", length(signature1)), rep("FOXP3", length(signature2))), stringsAsFactors = FALSE)
# Perform PCA
res=PCA(mtx_selected, quali.sup = ncol(mtx_selected), graph=T)
plot(res)
