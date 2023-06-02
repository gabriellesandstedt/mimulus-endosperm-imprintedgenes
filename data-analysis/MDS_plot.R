# Set working directory
setwd("~/Dropbox/research_UGA/working_rnaseq_folder")

# Load required packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("edgeR")
library(edgeR)

# Read data from file
data <- read.delim("Mopen_IM62v3_HTSeq_gene_counts_2023-05-11_sorted.txt", sep = "\t")

# Print first few rows and dimensions of data
head(data)
dim(data)

# Create edgeR object with counts and gene information
y <- DGEList(counts = data[, 2:19], genes = data[, 1:1])

# Perform normalization steps independent of experimental design
y <- calcNormFactors(y)

# Print scaling factors, which should be similar among all samples
y$samples

# Remove low expressed genes that do not have at least 1 CPM in at least 3 samples
isexpr <- rowSums(cpm(y) > 1) >= 3
y <- y[isexpr, , keep.lib.sizes = FALSE]

# Print resulting object
y
plotMDS(y)

# Custom MDS plot
points <- c(19, 19, 15, 3, 15, 3, 15, 3, 15, 3, 19, 19, 19, 19, 4, 15, 19, 4, 15, 3, 19, 15, 3)
colors <- rep(c("dodgerblue4", "dodgerblue", "darkorange", "darkorange3", "darkorange", "orange", "darkorange3", "orange", "dodgerblue", "dodgerblue4", "blue4", "darkorange3", "dodgerblue", "blue", "darkorange", "orange", "dodgerblue4", "darkorange3", "orange3"), 2)
plotMDS(y, pch = points, col = colors, cex = 3.3)
