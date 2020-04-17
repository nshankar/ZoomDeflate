# Set working directory, load libraries, and get source.
setwd('~/Documents/Projects/ZoomDeflate/')
library(ggplot2)
library(dplyr)
library(hrbrthemes)
source('ALRA/alra.R')
source('ALRA/jp_utilities.R')
set.seed(41) # Some ML projects make this seed a hyper-parameter


# Currently cheating by using the true dropouts as a mask
# In the future, use predicted dropouts as a mask
### PATHNAME is currently for Jeremy's computer, fix soon ###
PATHNAME = "~/Documents/Projects/ZoomDeflate/SplatGenData/5_groups_1000_cells_5000_genes/"
MASK_PATHNAME = paste(PATHNAME,"dropouts.csv", sep="")
DATA_PATHNAME = paste(PATHNAME,"counts.csv", sep="")
TRUE_PATHNAME = paste(PATHNAME,"true_counts.csv",sep="")

# load matrices in format: rows are genes, columns are cells
mask <- as.matrix(read.csv(MASK_PATHNAME, header=FALSE, sep=" "))
data <- as.matrix(read.csv(DATA_PATHNAME, header=FALSE, sep=" "))
truth <- as.matrix(read.csv(TRUE_PATHNAME, header=FALSE, sep=" "))

# take transpose b/c that's what ALRA likes 
mask <- t(mask)
data <- t(data)
truth <- t(truth)

# turn mask into logical array for logical indexing 
mask <- mask == 1

# Collect row sums to invert the normalize_data operation
row_sums_data = myRowSums(data)

# Library and log normalize the data
A_norm <- normalize_data(data)

# Choose k (# of singular values in the approximation) by measuring when they get smol. 
k_choice <- choose_k(A_norm)
# print(k_choice$k)

# complete matrix using ALRA
A_norm_completed <- alra(A_norm,k=k_choice$k)[[3]]
#A_norm_completed <- alra(A_norm,k=5)[[3]]

# invert the normalize_data operation
A_norm_completed <-  unnormalize_rows(A = A_norm_completed, og_row_sums = row_sums_data)


zero_stats <- zero_quality_stats(mask, truth, recon=A_norm_completed)
print(zero_stats)

RMSE_stats <- RMSE_for_sc(mask, truth, data, recon= A_norm_completed)  
print(RMSE_stats)



# # Calculate biological and technical zeros
# bio_zeros_mask <- (truth == 0)
# tech_zeros_mask <- (truth != 0) & (mask)
# 
# frac_bio_zeros <- sum(bio_zeros_mask) / length(truth)
# frac_tech_zeros <- sum(tech_zeros_mask) / length(truth)

# # Get some values for plotting
# reconned_dropouts <- A_norm_completed[tech_zeros_mask]
# truth_dropouts <- truth[tech_zeros_mask]
# reconned_nondropped <- A_norm_completed[!tech_zeros_mask]
# truth_nondropped <- truth[!tech_zeros_mask]
#

# 
# # histogram of dropout entries
# data <- data.frame(
#   type = c( rep("Reconned dropouts", length(reconned_dropouts)), rep("True dropouts", length(truth_dropouts))),
#   value = c( reconned_dropouts, truth_dropouts )
# )
# 
# p1 <- data %>%
#   ggplot( aes(x=value, fill=type)) +
#   geom_histogram( color="#e9ecef", alpha=0.6, position = 'identity',bins=100) +
#   scale_fill_manual(values=c("#69b3a2", "#404080")) +
#   theme_ipsum() +
#   labs(fill="") 
# 
# # histogram of nondropout entries
# data2 <- data.frame(
#   type = c( rep("Reconned nondropped", length(reconned_nondropped)), rep("True nondropped", length(truth_nondropped))),
#   value = c( reconned_nondropped, truth_nondropped )
# )
# 
# p2 <- data2 %>%
#   ggplot( aes(x=value, fill=type)) +
#   geom_histogram( color="#e9ecef", alpha=0.6, position = 'identity',bins=100) +
#   scale_fill_manual(values=c("#69b3a2", "#404080")) +
#   theme_ipsum() +
#   labs(fill="") 
# 
# # histogram of all entries
# data3 <- data.frame(
#   type = c( rep("Truth", length(truth)), rep("Reconned", length(A_norm_completed) )),
#   value = c( truth, A_norm_completed) )
# 
# p3 <- data3 %>%
#   ggplot( aes(x=value, fill=type)) +
#   geom_histogram( color="#e9ecef", alpha=0.6, position = 'identity',bins=100) +
#   scale_fill_manual(values=c("#69b3a2", "#404080")) +
#   theme_ipsum() +
#   labs(fill="") 
# 
# print("I'm done running")