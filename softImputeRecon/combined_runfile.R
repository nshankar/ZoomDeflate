# Nikhil's wd
setwd('/Users/nikhil/Documents/College/Math 651/ZoomDeflate/')

# Jeremy's wd
# setwd('~/Documents/Projects/ZoomDeflate/')

# Load some libraries 
library(softImpute)
library(ggplot2)
library(dplyr)
#library(hrbrthemes)
source('ALRA/alra.R')
source('ALRA/jp_utilities.R')

set.seed(43) # Some ML projects make this seed a hyper-parameter

# Data sets to run code upon
nGroups <- c(5)  #c(2, 5, 10)
nCells = c(1000) #c(1000, 10000)
nGenes = c(5000) #c(5000, 1000)
# Store RMSE results
softImpute_dict <- hash()
alra_dict <- hash()

# Loop structure depends heavily on directory structure
for (size in nGroups) {
  for (index in 1:length(nCells)){
    PATHNAME = paste("SplatGenData/", size , "_groups_", nCells[index], "_cells_", nGenes[index], "_genes/", sep="")
    MASK_PATHNAME = paste(PATHNAME, "dropouts.csv", sep="")
    DATA_PATHNAME = paste(PATHNAME, "counts.csv", sep="")
    TRUE_PATHNAME = paste(PATHNAME, "true_counts.csv", sep="")
    ID = paste(size, nCells[index], nGenes[index]) # for hashing
    
    # load matrices in format: rows are genes, columns are cells
    mask <- as.logical(as.matrix(read.csv(MASK_PATHNAME, header=FALSE, sep=" "), nrow = nCells[index], ncol = nGenes[index]))
    data <- as.matrix(read.csv(DATA_PATHNAME, header=FALSE, sep=" "), nrow = nCells[index], ncol = nGenes[index])
    truth <- as.matrix(read.csv(TRUE_PATHNAME, header=FALSE, sep=" "), nrow = nCells[index], ncol = nGenes[index])
    
    # some general preprocessing
    all_zeros_mask <- (data == 0)
    col_sums_data <- myColSums(data)
    data.normalized <- t(normalize_data(t(data)))
    
    ########### Run softImpute ###########
    # softImpute parameters
    lambda_= 1.9  # regularization weight emphasizing the nuclear norm
    type_ = "als" # either als or svd, als is faster
    rank.max_ = 0 # not a necessary argument, 
                  # but an estimate could be useful
    
    # set zeros to NA (unobserved entries) for the softImpute algo
    data.normalized[all_zeros_mask] <- NA
    
    if (rank.max_) {
      fits = softImpute(data.normalized, rank.max = rank.max_,
                        lambda = lambda_, trace=TRUE, type=type_)
    } else {
      fits = softImpute(data.normalized, lambda = lambda_,
                        trace=TRUE, type=type_)
    }
    
    # compute normalized output, then rescale so it is comparable with true_counts
    output_sI <- complete(data.normalized, fits)
    output_sI <-  unnormalize_data(output_sI, col_sums_data)
    
    # set zeros back to 0 in data (so that alra can run)
    data.normalized[all_zeros_mask] <- 0
    
    # compute some statistics
    zero_stats_sI <- zero_quality_stats(mask, truth, recon=output_sI)
    print(zero_stats_sI)
    
    RMSE_stats_sI <- RMSE_for_sc(mask, truth, data, recon=output_sI)
    softImpute_dict[[ID]] <- RMSE_stats_sI
    print(RMSE_stats_sI)
    
    
    ########### RUN ALRA ###########
    # take transpose b/c that's what ALRA enjoys 
    mask <- t(mask)
    A_norm <- t(data.normalized) #rename normalized data
    truth <- t(truth)
    row_sums_data <- col_sums_data
    
    # Choose k (# of singular values in the approximation) by measuring when they get smol. 
    k_choice <- choose_k(A_norm)
    
    # complete matrix using ALRA
    output_ALRA <- alra(A_norm,k=k_choice$k)[[3]]
    
    # invert the normalize_data operation
    output_ALRA <-  unnormalize_rows(A = output_ALRA, og_row_sums = row_sums_data)
    
    # compute some statistics
    zero_stats_ALRA <- zero_quality_stats(mask, truth, recon=output_ALRA)
    print(zero_stats_ALRA)
    
    RMSE_stats_ALRA <- RMSE_for_sc(mask, truth, data, recon= output_ALRA)
    alra_dict[[ID]] <- RMSE_stats_ALRA
    print(RMSE_stats_ALRA)
  }
}
  
  #### RMSE-hacking
  # 