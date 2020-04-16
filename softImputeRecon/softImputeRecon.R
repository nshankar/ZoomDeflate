library(softImpute)
source("/Users/nikhil/Documents/College/Math 651/ZoomDeflate/ALRA/alra.R")



unnormalize_data <- function(A, og_col_sums) {
  # invert the function normalize_data from alra.R
  A <- (exp(A) - 1)/1E4
  A <- sweep(A, 2, col_sums_data, '*')
  return(A)
}

myColSums <- function(A) {
  # get column sums, meant to work specifically with normalize_data from alra.R
  col_sums= colSums(A)
  if (any(col_sums == 0)) {
    toRemove <- which(col_sums == 0)
    data <- data[,-toRemove]
    col_sums <- col_sums[-toRemove]
  }
  return(col_sums)
}


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


mask <- as.matrix(read.table(MASK_PATHNAME), nrow = nCells, ncol = nGenes)
mask <- as.logical(mask)
true_counts <- as.matrix(read.table(DATA_PATHNAME), nrow = nCells, ncol = nGenes)

#normalize in the correct order
data <- true_counts
data[mask] = 0
# Collect column sums to invert the normalize_data operation
col_sums_data = myColSums(data)
data <- t(normalize_data(t(data)))
data[mask] = NA


if (rank.max_) {
  fits = softImpute(data, rank.max = rank.max_,
                    lambda = lambda_, trace=TRUE, type=type_)
} else {
  fits = softImpute(data, lambda = lambda_,
                    trace=TRUE, type=type_)
}

# compute normalized output, then rescale so it is comparable with true_counts
output <- complete(data, fits)
output <-  unnormalize_data(output, col_sums_data)
