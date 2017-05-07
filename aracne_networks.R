# https://bioconductor.org/packages/devel/data/experiment/vignettes/aracne.networks/inst/doc/aracne.networks.pdf
options(stringsAsFactors = FALSE)
library(aracne.networks)
library(dplyr)
library(networkD3)
library(igraph)

# All available cancers
data(package="aracne.networks")$results[, "Item"]

# Full network, top 10
data(regulonbrca)
tmp <- write.regulon(regulonblca, n = 10) %>% data.frame

# Network for a specific regulator. Entrez IDs
write.regulon(regulonbrca, regulator="6382", file = "aracne.txt") # SDC1 - 6382

gene <- "4605"
networkData <- data.frame(regulonbrca$`4605`)
networkData <- data.frame(Source = gene, Target = rownames(networkData), MoA = networkData$tfmode, likelihood = networkData$likelihood)
simpleNetwork(networkData, zoom = TRUE)

# Attempt to do directed graph. If MoA is negative, swap source and target
for (i in 1:nrow(networkData)) {
  if (networkData[i, "MoA"] < 0) {
    networkData[i, 1:2] <- as.character(networkData[i, 2:1])
  }
}
# Order by source
networkData <- networkData[order(networkData$Source), ]
# Create notes dataframe
networkNodes <- data.frame(name = unique(c(networkData$Source, networkData$Target), group = 1))
networkNodes$group[ networkNodes$name == gene ] <- 2 # Set group for the gene of interest

# Create network links, with value as an inverse of the likelihood
networkLinks <- data.frame(source = networkData$Source, target = networkData$Target, value = 1 / networkData$likelihood)

# Doesn't work - need links as indices, not IDs
# forceNetwork(Links = networkLinks, Nodes = networkNodes, Source = "source",
#              Target = "target", Value = "value", NodeID = "name",
#              Group = "group", opacity = 0.4, zoom = TRUE)

# Using igraph workaround
net   <- graph_from_data_frame(d = networkLinks, directed = TRUE, vertices = networkNodes)
net3D <- igraph_to_networkD3(net, group = networkNodes$group)

forceNetwork(Links = net3D$links, Nodes = net3D$nodes, Source = "source",
             Target = "target", NodeID = "name", Group = "group", Value = "value")

load(file = "/Users/mdozmorov/Documents/Work/VCU_work/Mark/disease-coherence/data/9606.protein.links.v10.rda")
load(file = "/Users/mdozmorov/Documents/Work/VCU_work/Mark/disease-coherence/data/biogrid.rda")
load(file = "/Users/mdozmorov/Documents/Work/VCU_work/Mark/disease-coherence/data/i2d.2_9.Public.HUMAN.tab.rda")
head(my_adj_list)
selected <- my_adj_list[ (my_adj_list$from %in% "SDC1") | (my_adj_list$to %in% "SDC1"), , drop = FALSE ]
dim(selected)
head(selected)
selected <- selected[complete.cases(selected), ]
selected$to[ grep("CXCL8", selected$to) ] %>% unique %>% sort
selected$from[ grep("CXCL8", selected$from) ] %>% unique %>% sort
simpleNetwork(selected, Source = "from", Target = "to", fontSize = 10, linkColour = "black", nodeColour = "blue", opacity = 0.8, zoom = TRUE)


net <- graph_from_data_frame(d = selected, directed = TRUE)
plot(net, edge.arrow.size=.6, vertex.label.family = "Arial")
