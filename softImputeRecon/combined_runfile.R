# Nikhil's wd
setwd('/Users/nikhil/Documents/College/Math 651/ZoomDeflate/')

# Jeremy's wd
# setwd('~/Documents/Projects/ZoomDeflate/')

# Load some libraries 
library(softImpute)
library(ggplot2)
library(dplyr)
library(hash)
library(Rtsne)
library(plotly)
#library(hrbrthemes)
source('ALRA/alra.R')
source('ALRA/jp_utilities.R')

set.seed(43) # Some ML projects make this seed a hyper-parameter

# Data sets to run code upon
nGroups <- c(2, 5, 10)
nCells <- c(1000, 10000)
nGenes <- c(5000, 1000)

# Store RMSE results
softImpute_dict <- hash()
alra_dict <- hash()
merged_dict <- hash()
softImpute_thresh_dict <- hash()
alra_RMSE_hack_dict <- hash()

alra_ranks <- hash()

observed_dict <- hash()

#Should I do t-SNE
switch_tSNE = FALSE
switch_writeOutput = FALSE

# Loop structure depends heavily on directory structure
for (size in nGroups) {
  for (index in 1:length(nCells)){
    PATHNAME = paste("SplatGenData/", size , "_groups_", nCells[index], "_cells_", nGenes[index], "_genes/", sep="")
    MASK_PATHNAME = paste(PATHNAME, "dropouts.csv", sep="")
    DATA_PATHNAME = paste(PATHNAME, "counts.csv", sep="")
    TRUE_PATHNAME = paste(PATHNAME, "true_counts.csv", sep="")
    ID = paste("(", size, ", ", nCells[index], ", ", nGenes[index], ")", sep="") # for hashing
    
    # load matrices in format: rows are genes, columns are cells
    mask <- as.matrix(read.csv(MASK_PATHNAME, header=FALSE, sep=" ")) # , nrow = nCells[index], ncol = nGenes[index])
    data <- as.matrix(read.csv(DATA_PATHNAME, header=FALSE, sep=" "), nrow = nCells[index], ncol = nGenes[index])
    truth <- as.matrix(read.csv(TRUE_PATHNAME, header=FALSE, sep=" "), nrow = nCells[index], ncol = nGenes[index])
    mask <- (mask == 1)
    
    zero_stats_observed <- zero_quality_stats(mask, truth, recon=data)
    RMSE_stats_observed <- RMSE_for_sc(mask, truth, data, data)
    observed_dict[[ID]] <- c(RMSE_stats_observed, zero_stats_observed)
    
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
    output_sI <-  round(unnormalize_data(output_sI, col_sums_data))
    
    # set zeros back to 0 in data (so that alra can run)
    data.normalized[all_zeros_mask] <- 0
    
    # compute some statistics
    zero_stats_sI <- zero_quality_stats(mask, truth, recon=output_sI)
    print(zero_stats_sI)
    
    RMSE_stats_sI <- RMSE_for_sc(mask, truth, data, recon=output_sI)
    softImpute_dict[[ID]] <- c(RMSE_stats_sI,zero_stats_sI)
    #print(RMSE_stats_sI)
    
    if (switch_writeOutput) {
      write.csv(output_sI,paste(PATHNAME, "output_sI.csv",sep=""),row.names=FALSE)
    }
    if (switch_tSNE) {
      tSNE_sI <- tSNE_2(output_sI)
      write.csv(tSNE_sI,paste(PATHNAME, "tSNE_sI.csv",sep=""),row.names=FALSE)
    }
    
    ########### RUN ALRA ###########
    # take transpose b/c that's what ALRA enjoys 
    A_norm <- t(data.normalized) #rename normalized data
    row_sums_data <- col_sums_data
    
    # Choose k (# of singular values in the approximation) by measuring when they get smol. 
    k_choice <- choose_k(A_norm)
    
    alra_ranks[[ID]] <- k_choice$k
    
    # complete matrix using ALRA
    output_ALRA <- alra(A_norm,k=k_choice$k)[[3]]
    
    # invert the normalize_data operation
    output_ALRA <-  round(unnormalize_rows(A = output_ALRA, og_row_sums = row_sums_data))
    
    # take transpose for better handling
    output_ALRA <- t(output_ALRA)
    
    # compute some statistics
    zero_stats_ALRA <- zero_quality_stats(mask, truth, recon=output_ALRA)
    print(zero_stats_ALRA)
    
    RMSE_stats_ALRA <- RMSE_for_sc(mask, truth, data, recon= output_ALRA)
    alra_dict[[ID]] <- c(RMSE_stats_ALRA,zero_stats_ALRA)

    print(RMSE_stats_ALRA)
    if (switch_writeOutput) {
      write.csv(output_ALRA,paste(PATHNAME, "output_ALRA.csv",sep=""),row.names=FALSE)
    }
    if (switch_tSNE) {
      tSNE_ALRA <- tSNE_2(output_ALRA)
      write.csv(tSNE_ALRA,paste(PATHNAME, "tSNE_ALRA.csv",sep=""),row.names=FALSE)
    }
    ########### METHOD 3:  ###########
    # output_merged <- output_sI
    # output_merged[output_ALRA==0] <- 0
    # 
    # zero_stats_merged <- zero_quality_stats(mask, truth, recon=output_merged)
    # print(zero_stats_merged)
    # 
    # RMSE_stats_merged <- RMSE_for_sc(mask, truth, data, recon= output_merged)
    # print(RMSE_stats_merged)
    # 
    # merged_dict[[ID]] <- c(RMSE_stats_merged,zero_stats_merged)
    # 
    # if (switch_writeOutput) {
    #   write.csv(output_merged,paste(PATHNAME, "output_merged.csv",sep=""),row.names=FALSE)
    # }
    # if (switch_tSNE) {
    #   tSNE_merged <- tSNE_2(output_merged)
    #   write.csv(tSNE_merged,paste(PATHNAME, "tSNE_merged.csv",sep=""),row.names=FALSE)
    # }
    # ########### METHOD 4: softImpute thresholded ###########
    # output_sI_thresh <- output_sI
    # output_sI_thresh[output_sI < 0.5] <- 0
    # 
    # zero_stats_thresh <- zero_quality_stats(mask, truth, recon=output_sI_thresh)
    # print(zero_stats_thresh)
    # 
    # RMSE_stats_thresh <- RMSE_for_sc(mask, truth, data, recon= output_sI_thresh)
    # print(RMSE_stats_thresh)
    # 
    # softImpute_thresh_dict[[ID]] <- c(RMSE_stats_thresh,zero_stats_thresh)
    # 
    # write.csv(output_sI_thresh,paste(PATHNAME, "output_sI_thresh.csv",sep=""),row.names=FALSE)
    # 
    # tSNE_sI_thresh <- tSNE_2(output_sI_thresh)
    # write.csv(tSNE_sI_thresh,paste(PATHNAME, "tSNE_sI_thresh.csv",sep=""),row.names=FALSE)
    # 
    ########### METHOD 5: ALRA RMSE Hack ###########
    output_ALRA_hack <- output_ALRA
    output_ALRA_hack[!all_zeros_mask] <- data[!all_zeros_mask]
    
    # compute some statistics
    zero_stats_ALRA_hack <- zero_quality_stats(mask, truth, recon=output_ALRA_hack)
    print(zero_stats_ALRA_hack)
    
    RMSE_stats_ALRA_hack <- RMSE_for_sc(mask, truth, data, recon= output_ALRA_hack)
    print(RMSE_stats_ALRA_hack)
    alra_RMSE_hack_dict[[ID]] <- c(RMSE_stats_ALRA_hack,zero_stats_ALRA_hack)
    
    if(switch_writeOutput) {
      write.csv(output_ALRA_hack,paste(PATHNAME, "output_ALRA_hack.csv",sep=""),row.names=FALSE)
    }
    if(switch_tSNE) {
      tSNE_ALRA_hack <- tSNE_2(output_ALRA_hack)
      write.csv(tSNE_ALRA_hack,paste(PATHNAME, "tSNE_ALRA_hack.csv",sep=""),row.names=FALSE)
    }
  }
}
#   
# 
# output_sI_hack <- output_sI 
# output_sI_hack[!all_zeros_mask] <- data[!all_zeros_mask]
# 
# # # compute some statistics
# zero_stats_sI_hack <- zero_quality_stats(mask, truth, recon=output_sI_hack)
# print(zero_stats_sI_hack)
# 
# RMSE_stats_sI_hack <- RMSE_for_sc(mask, truth, data, recon= output_sI_hack)
# print(RMSE_stats_sI_hack)
# 
# output_ALRA_hack <- output_ALRA
# output_ALRA_hack[!all_zeros_mask] <- data[!all_zeros_mask]
# 
# # compute some statistics
# zero_stats_ALRA_hack <- zero_quality_stats(mask, truth, recon=output_ALRA_hack)
# print(zero_stats_ALRA_hack)
# 
# RMSE_stats_ALRA_hack <- RMSE_for_sc(mask, truth, data, recon= output_ALRA_hack)
# print(RMSE_stats_ALRA_hack)


myZeroQualBarplots()

##################### \begin{RMSE Chart} ####################################
library(reactable)

RMSE_dicts <- c(softImpute_dict, alra_RMSE_hack_dict, alra_dict, observed_dict)
RMSE_matrix <- matrix(NA, nrow = 6, ncol = 8)

observed_RMSE

alra_hack_string <- "ALRA-Alt" #find a better name for this
row_names<-rep("", 6)

for (j in 1:length(nCells)){
  for (i in 1:length(nGroups)) {
    ID = paste("(", nGroups[i], ", ", nCells[j], ", ", nGenes[j], ")", sep="")
    row <- i + (j-1)*length(nGroups)
    col <- 1
    row_names[row] <- ID
    for (dict in RMSE_dicts) {
      RMSE_matrix[row, col] <- dict[[ID]]$RMSE_all
      col <- col + 1
      RMSE_matrix[row, col] <- dict[[ID]]$RMSE_dropouts
      col <- col + 1
    }
  }
}


RMSE_matrix <- round(RMSE_matrix, digits=2) # for readability
RMSE_data <- as.data.frame(RMSE_matrix, row.names = row_names)

### Trying some 538 stuff ###
maxWidth_ = 75
# Ce n'est pas une pal orange
orange_pal <- function(x) rgb(colorRamp(c("ivory", "lightcoral"))(x), maxColorValue = 255)

mystyle <- function(value) {
  normalized <- sqrt((value - min(RMSE_matrix)) / (max(RMSE_matrix) - min(RMSE_matrix)))
  color <- orange_pal(normalized)
  list(background = color)
}

reactable(
  RMSE_data, rownames=TRUE,
  style = list(fontFamily = "Work Sans, sans-serif", fontSize = "14px"),
  columns = list(
    .rownames = colDef(name = "(# Groups,    # Cells, # Genes)",
                       headerStyle = list(fontSize = "12px")),
    V1 = colDef(name = "All", style = mystyle, 
                maxWidth = maxWidth_, headerStyle = list(fontSize = "12px")),
    V2 = colDef(name = "Dropouts", style = mystyle, 
                maxWidth = maxWidth_, headerStyle = list(fontSize = "12px")),
    V3 = colDef(name = "All", style = mystyle, 
                maxWidth = maxWidth_, headerStyle = list(fontSize = "12px")),
    V4 = colDef(name = "Dropouts", style = mystyle, 
                maxWidth = maxWidth_, headerStyle = list(fontSize = "12px")),
    V5 = colDef(name = "All", style = mystyle, 
                maxWidth = maxWidth_, headerStyle = list(fontSize = "12px")),
    V6 = colDef(name = "Dropouts", style = mystyle, 
                maxWidth = maxWidth_, headerStyle = list(fontSize = "12px")),
    V7 = colDef(name = "All", style = mystyle, 
                maxWidth = maxWidth_, headerStyle = list(fontSize = "12px")),
    V8 = colDef(name = "Dropouts", style = mystyle, 
                maxWidth = maxWidth_, headerStyle = list(fontSize = "12px"))
  ),
  columnGroups = list(
    colGroup(name = "Data Sets", columns = c(".rownames"), align = "left"),
    colGroup(name = "ZoomDeflate", columns = c("V1", "V2"), align="right"),
    colGroup(name = alra_hack_string, columns = c("V3", "V4"), align="right"),
    colGroup(name = "ALRA", columns = c("V5", "V6"), align="right"),
    colGroup(name="Observed", columns=c("V7", "V8"), align="right")
  )
)
################## \end{RMSE Chart} #####################################

# # tSNE accepts objects as rows, dimensions as columns 
# # I don't know what the normalization is. 
# output_sI_normed <- normalize_input(t(output_sI_hack)) 
# output_ALRA_normed <- normalize_input(t(output_ALRA_hack))
# output_merged_normed <- normalize_input(t(output_merged))
# 
# # performs tSNE + PCA 
# low_dim_rep_sI <- Rtsne(output_sI_normed)$Y
# low_dim_rep_ALRA <- Rtsne(output_ALRA_normed)$Y
# low_dim_rep_merged <- Rtsne(output_merged_normed)$Y
# 
# 
# sI_fr <- as.data.frame(low_dim_rep_sI)
# ALRA_fr <- as.data.frame(low_dim_rep_ALRA)
# merged_fr <- as.data.frame(low_dim_rep_merged)
# 
# 
# fig1 <- plot_ly(data = sI_fr,x=~V1,y=~V2,mode='markers')
# fig2 <- plot_ly(data = ALRA_fr,x=~V1,y=~V2,mode='markers')
# fig3 <- plot_ly(data = merged_fr,x=~V1,y=~V2,mode='markers')
# fig1
# fig2
# fig3


####
# write the matrices to csv
# load them into MATLAB/Python and do a clustering analysis
# (or do that natively in R)
# with goal of computing consistency, ARI, something like this.
# Look at tSNE plot
#
# 

#### RMSE-hacking
# data <- t(data)
# all_zeros_mask <- t(all_zeros_mask)
# output_ALRA[!all_zeros_mask] <- data[!all_zeros_mask]
# data <- t(data)
# all_zeros_mask <- t(all_zeros_mask)
# 
# # compute some statistics
# zero_stats_ALRA2 <- zero_quality_stats(t(mask), t(truth), recon=output_ALRA)
# print(zero_stats_ALRA2)
# 
# RMSE_stats_ALRA2 <- RMSE_for_sc(t(mask), t(truth), t(data), recon= output_ALRA)  
# print(RMSE_stats_ALRA2)
# 
# output_sI[!all_zeros_mask] <- data[!all_zeros_mask]
# 
# # compute some statistics
# zero_stats_sI2 <- zero_quality_stats(mask, truth, recon=output_sI)
# print(zero_stats_sI2)
# 
# RMSE_stats_sI2 <- RMSE_for_sc(mask, truth, data, recon= output_sI)  
# print(RMSE_stats_sI2)
