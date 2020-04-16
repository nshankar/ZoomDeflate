# Nikhil's wd
setwd('/Users/nikhil/Documents/College/Math 651/ZoomDeflate/softImputeRecon')

# Jeremy's wd
#setwd('~/Documents/Projects/ZoomDeflate/softImputeRecon')

library(softImpute)
library(hash)
source("../ALRA/alra.R")
source("../ALRA/jp_utilities.R")

# Currently cheating by using the true dropouts as a mask
# In the future, use predicted dropouts as a mask

# These inputs depend heavily on the directory structure
nGroups <- c(2, 5, 10)
nCells = c(1000, 10000)
nGenes = c(5000, 1000)

stats_dict <- hash()

for (size in nGroups) {
  for (index in 1:length(nCells)){
    PATHNAME = paste("../SplatGenData/", size , "_groups_", nCells[index], "_cells_", nGenes[index], "_genes/", sep="")
    MASK_PATHNAME = paste(PATHNAME, "dropouts.csv", sep="")
    DATA_PATHNAME = paste(PATHNAME, "true_counts.csv", sep="")
    
    
    set.seed(1011)# Some ML projects make this seed a hyper-parameter
    lambda_= 1.9  # regularization weight emphasizing the nuclear norm
    type_ = "als" # either als or svd, als is faster
    rank.max_ = 0 # not a necessary argument, 
                  # but an estimate could be useful
    
    
    mask <- as.matrix(read.table(MASK_PATHNAME), nrow = nCells[index], ncol = nGenes[index])
    mask <- as.logical(mask)
    true_counts <- as.matrix(read.table(DATA_PATHNAME), nrow = nCells[index], ncol = nGenes[index])
    
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
    
    zero_stats <- zero_quality_stats(mask, truth=true_counts, recon=output)
    print(zero_stats)
    
    data[mask] <- 0
    data <- unnormalize_data(data, col_sums_data)
    RMSE_stats <- RMSE_for_sc(mask, truth=true_counts, data, output)
    id <- paste(size, nCells[index], nGenes[index])
    stats_dict[[id]] <- RMSE_stats
    #print(RMSE_stats)
  }
}
