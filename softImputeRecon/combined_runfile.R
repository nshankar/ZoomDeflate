# Nikhil's wd
# setwd('/Users/nikhil/Documents/College/Math 651/ZoomDeflate/softImputeRecon/')

# Jeremy's wd
setwd('~/Documents/Projects/ZoomDeflate/')

# Load some libraries 
library(softImpute)
library(ggplot2)
library(dplyr)
library(hrbrthemes)
source('ALRA/alra.R')
source('ALRA/jp_utilities.R')

# Filenames of input 
data_path = "SplatGenData/5_groups_1000_cells_5000_genes/"
MASK_PATHNAME = paste(data_path,"dropouts.csv", sep="")
DATA_PATHNAME = paste(data_path,"counts.csv", sep="")
TRUE_PATHNAME = paste(data_path,"true_counts.csv",sep="")

# load matrices in format: rows are genes, columns are cells
mask <- as.matrix(read.csv(MASK_PATHNAME, header=FALSE, sep=" "))
data <- as.matrix(read.csv(DATA_PATHNAME, header=FALSE, sep=" "))
truth <- as.matrix(read.csv(TRUE_PATHNAME, header=FALSE, sep=" "))

# turn mask into logical array for logical indexing 
mask <- mask == 1
# get mask of all zeros
all_zeros_mask <- data == 0

# Parameters
nCells = NCOL(mask)
nGenes = NROW(mask)

set.seed(41) # Some ML projects make this seed a hyper-parameter
lambda_= 1.9  # regularization weight emphasizing the nuclear norm
type_ = "als" # either als or svd, als is faster
rank.max_ = 0 # not a necessary argument, 
              # but an estimate could be useful

### RUN ALRA ### 
# take transpose b/c that's what ALRA takes 
mask <- t(mask)
data <- t(data)
truth <- t(truth)

# Collect row sums to invert the normalize_data operation
row_sums_data = myRowSums(data)
data <- normalize_data(data)

# Library/log normalize the data
A_norm <- normalize_data(data)

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
print(RMSE_stats_ALRA)

### RUN SOFTIMPUTE ### 
# Undo the transpose b/c that's what softImpute takes 
mask <- t(mask)
data <- t(data)
truth <- t(truth) 

# note that data is already scaled appropriately. 
# set zeros to NA (unobserved entries)
data[all_zeros_mask] <- NA

# collect column sums (row sums of transpose)
col_sums_data <- row_sums_data

#  Run soft-impute 
if (rank.max_) {
  fits = softImpute(data, rank.max = rank.max_,
                    lambda = lambda_, trace=TRUE, type=type_)
} else {
  fits = softImpute(data, lambda = lambda_,
                    trace=TRUE, type=type_)
}

# compute normalized output, then rescale so it is comparable with true_counts
output_sI <- complete(data, fits)
output_sI <-  unnormalize_data(output_sI, col_sums_data)

# set zeros back to 0 in data
data[all_zeros_mask] <- 0
data <- unnormalize_data(data, col_sums_data)

# compute some statistics
zero_stats_sI <- zero_quality_stats(mask, truth, recon=output_sI)
print(zero_stats_sI)

RMSE_stats_sI <- RMSE_for_sc(mask, truth, data, recon= output_sI)  
print(RMSE_stats_sI)


#### RMSE-hacking
# 