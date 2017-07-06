set.seed(1)
# A vector of survivalCpG and unsurvivalCpG sites
y <- c(rep("survivalCpG", 20), rep("unsurvivalCpG", 20))
y
# Vectors of overlap assignment
WithinGene <- Distant <- rep(0, 40) # Initialize with 0 - no overlaps
WithinGene[rbinom(60, 40, 0.4)]   <- 1 # Make most "survivalCpG" regions overlap "WithinGene" annotation
Distant[rbinom(60, 40, 0.6)] <- 1 # Make most "unsurvivalCpG" regions overlap "Distant" annotation
# Assemble all in a data frame
mtx <- data.frame(y = y, WithinGene = WithinGene, Distant = Distant, stringsAsFactors = FALSE)
mtx

# Get pieces of contingency table
annotation <- mtx$WithinGene # Select which annotation to test
# "Y" or "N" suffices correspond to "Yes" or "No"
survivalCpGY_annotationY <- sum(y == "survivalCpG"   & annotation == 1) # How many are "survivalCpG" and "1" (overlap annotation)
survivalCpGY_annotationN <- sum(y == "survivalCpG"   & annotation == 0)
survivalCpGN_annotationY <- sum(y == "unsurvivalCpG" & annotation == 1)
survivalCpGN_annotationN <- sum(y == "unsurvivalCpG" & annotation == 0)
# Assemble contingency table
mtx_2x2 <- matrix(c(survivalCpGY_annotationY, survivalCpGY_annotationN,
                    survivalCpGN_annotationY, survivalCpGN_annotationN),
                  nrow = 2, dimnames = list(c("annotationY", "annotationN"), c("survivalCpGY", "survivalCpGN")))
mtx_2x2
# Fisher's test whether survivalCpG regions are more significantly overlap with the selected annotation
fisher.test(mtx_2x2)
