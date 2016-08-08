# https://bioconductor.org/packages/devel/data/experiment/vignettes/aracne.networks/inst/doc/aracne.networks.pdf

library(aracne.networks)
data(package="aracne.networks")$results[, "Item"]

data(regulonblca)
write.regulon(regulonblca, n = 10)

data(regulonblca)
write.regulon(regulonblca, regulator="399")
