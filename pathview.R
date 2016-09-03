library(pathview)
library(openxlsx)
degs <- read.xlsx("PAAD_CD31_hi51_lo49_DEGs.xlsx", cols = c(1, 2))
degs.genes <- degs$logFC
names(degs.genes) <- degs$Gene

pv.out <- pathview(gene.data = degs.genes, pathway.id = "04514", species = "hsa", gene.idtype = "SYMBOL", gene.annotpkg = "org.Hs.eg.db", out.suffix = "CD31")
