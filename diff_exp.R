# Prerequisites, prepared by the survival.R script
# expr - expression matrix separated by the high/low expression of the selected genes
# group - labeling of samples having high/low expression of the selected genes

# Reshape expression matrix
expr <- (t(expr))
colnames(expr) <- expr[1, ]
expr <- expr[-1, ]
class(expr) <- "numeric"
expr <- log2(expr)
boxplot(expr)

library(limma)
library(openxlsx)
library(MDmisc)
library(org.Hs.eg.db)

# Limma
design <- model.matrix(~0 + factor(group))
colnames(design) <- c("up", "lo")
fit <- lmFit(expr, design)
contrast.matrix <- makeContrasts(up-lo, levels = design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
topTable(fit2, coef=1, adjust="BH")

degs <- topTable(fit2, coef = 1, number = Inf, p.value = 0.01)

# A wrapper function to perform all functional enrichment analyses.
# Helper function to save non-empty results
save_res <- function(res, fileName = fileName, wb = wb, sheetName = "KEGG") {
  if (nrow(res) > 0) {
    openxlsx::addWorksheet(wb = wb, sheetName = sheetName)
    openxlsx::writeData(wb, res, sheet = sheetName)
    openxlsx::saveWorkbook(wb, fileName, overwrite = TRUE)
  }
}
# Create (or, load)  Excel file
fileName <- "fileName.xlsx"
wb <- openxlsx::createWorkbook(fileName) # openxlsx::loadWorkbook(fileName)
openxlsx::addWorksheet(wb = wb, sheetName = "DEGs")
openxlsx::writeData(wb, degs, sheet = "DEGs")
openxlsx::saveWorkbook(wb, fileName, overwrite = TRUE)

# Gene ontology, molecular function
res <- gene_enrichment(selected = rownames(degs), id="symbol", organism = "Hs", use="GO", ont="MF")
save_res(res, fileName, wb = wb, sheetName = "GOMF")
# Gene ontology, biological process 
res <- gene_enrichment(selected = rownames(degs), id="symbol", organism = "Hs", use="GO", ont="BP")
save_res(res, fileName, wb = wb, sheetName = "GOBP")
# Gene ontology, cellular component
res <- gene_enrichment(selected = rownames(degs), id="symbol", organism = "Hs", use="GO", ont="CC")
save_res(res, fileName, wb = wb, sheetName = "GOCC")
# KEGG canonical pathways
res <- gene_enrichment(selected = rownames(degs), id="symbol", organism = "Hs", use="KEGG")
save_res(res, fileName, wb = wb, sheetName = "KEGG")

# Plotting
index.to.plot <- order(group)
matrix.to.plot <- expr[rownames(degs)[1:50], index.to.plot]
genes.to.plot <- rownames(degs)[1:50]
group.to.plot <- group[index.to.plot]

NMF::aheatmap(matrix.to.plot, color=colorRampPalette(c('blue', 'gray', 'yellow'))(20), Colv = NA, Rowv = FALSE, hclust = "ward", scale = "row", annCol = group.to.plot, annColors = list(c("red", "blue")), labRow = genes.to.plot, fontsize = 10, cexRow = 10) # color="-RdYlBu"


