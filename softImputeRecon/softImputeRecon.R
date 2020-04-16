library(softImpute)
source("/Users/nikhil/Documents/College/Math 651/ZoomDeflate/ALRA/alra.R")

# Currently cheating by using the true dropouts as a mask
# In the future, use predicted dropouts as a mask
### PATHNAME is currently for Nikhil's computer, fix soon ###

nCells = 1000
nGenes = 5000

PATHNAME = "/Users/nikhil/Documents/College/Math 651/ZoomDeflate/SplatGenData/5_groups_1000_cells_5000_genes/"
MASK_PATHNAME = paste(PATHNAME, "dropouts.csv", sep="")
DATA_PATHNAME = paste(PATHNAME, "true_counts.csv", sep="")


set.seed(1011)# Some ML projects make this seed a hyper-parameter
lambda_= 1.9  # regularization weight emphasizing the nuclear norm
type_ = "als" # either als or svd, als is faster
rank.max_ = 0 # not a necessary argument, 
              # but an estimate could be useful

normalize_rows_by_noise = 0 #boolean


mask <- as.matrix(read.table(MASK_PATHNAME), nrow = nCells, ncol = nGenes)
mask <- as.logical(mask)
true_counts <- as.matrix(read.table(DATA_PATHNAME), nrow = nCells, ncol = nGenes)

#normalize in the correct order
data <- true_counts
true_counts <- t(normalize_data(t(true_counts)))
data[mask] = 0
data <- t(normalize_counts(t(data)))
data[mask] = NA


if (rank.max_) {
  fits = softImpute(data, rank.max = rank.max_,
                    lambda = lambda_, trace=TRUE, type=type_)
} else {
  fits = softImpute(data, lambda = lambda_,
                    trace=TRUE, type=type_)
}

output <- complete(data, fits)

diff_squared = (output - og_data)^2 #This isn't right.