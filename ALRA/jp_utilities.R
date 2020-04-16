# Inputs: mask (logical matrix with TRUE in locations of dropouts)
#         truth (matrix with true counts. Should have same scaling as recon)
#         recon (matrix with reconstructed counts)

# Outputs: fraction of true zeros preserved, fraction of fake zeros reconstructed
zero_quality_stats <- function(mask, truth, recon) {
  # Calculate biological and technical zeros
  bio_zeros_mask <- (truth == 0)
  tech_zeros_mask <- (truth != 0) & (mask)
  
  frac_bio_zeros <- sum(bio_zeros_mask) / length(truth)
  frac_tech_zeros <- sum(tech_zeros_mask) / length(truth)

  # Count biological zeros preserved
  bio_zeros_preserved <- sum(recon[bio_zeros_mask] == 0) 
  frac_bio_zeros_preserved <- bio_zeros_preserved / sum(bio_zeros_mask)

  # Count technical zeros reconstructed
  tech_zeros_reconned <- sum(recon[tech_zeros_mask] != 0) 
  frac_tech_zeros_reconned <- tech_zeros_reconned / sum(tech_zeros_mask)

return (list(frac_bio_zeros_preserved = frac_bio_zeros_preserved,
             frac_tech_zeros_reconned = frac_tech_zeros_reconned))
}


# Inputs: mask (logical matrix with TRUE in locations of dropouts)
#         truth (matrix with true counts. Should have same scaling as data/recon)
#         data (matrix with observed counts. should have same scaling as truth/recon)
#         recon (matrix with reconstructed counts.)
# Outputs: several RMSEs for different subsets of entries
RMSE_for_sc <- function(mask, truth, data, recon) {
  # Calculate biological and technical zeros
  tech_zeros_mask <- (truth != 0) & (mask)
  # Calculate RMSE for all values, 
  difference <- recon - truth
  RMSE_all <- sqrt(sum(difference^2) / (length(truth)))

  reconned_dropouts <- recon[tech_zeros_mask]
  truth_dropouts <- truth[tech_zeros_mask]
  RMSE_dropouts <- sqrt(sum((reconned_dropouts - truth_dropouts)^2) / sum(tech_zeros_mask))

  reconned_nondropped <- recon[!tech_zeros_mask]
  truth_nondropped <- truth[!tech_zeros_mask]
  RMSE_nondropped <- sqrt(sum((reconned_nondropped - truth_nondropped)^2) / sum(!tech_zeros_mask))
  
  zeros_mask <- data == 0
  reconned_zeros <- recon[zeros_mask]
  truth_zeros <- truth[zeros_mask]
  RMSE_zeros <- sqrt(sum((reconned_zeros-truth_zeros)^2) / sum(zeros_mask))
  
  reconned_nonz <- recon[!zeros_mask]
  truth_nonz <- truth[!zeros_mask]
  RMSE_nonz <- sqrt(sum((reconned_nonz - truth_nonz)^2) / sum(!zeros_mask))
  
  return (list(RMSE_all = RMSE_all,
             RMSE_dropouts = RMSE_dropouts,
             RMSE_nondropouts = RMSE_nondropped,
             RMSE_zeros = RMSE_zeros,
             RMSE_nonzeros = RMSE_nonz))
}
 
## ATT'n NIKHIL: should col_sums_data be og_col_sums in the sweep command? 
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

myRowSums <- function(A) {
  # get row sums, meant to work specifically with normalize_data from alra.R
  row_sums = rowSums(A)
  if (any(row_sums == 0)) {
    toRemove <- which(row_sums == 0)
    data <- data[-toRemove,]
    row_sums <- row_sums[-toRemove]
  }
  return(row_sums)
}

unnormalize_rows <- function(A, og_row_sums) {
  # invert the function normalize_data from alra.R
  A <- (exp(A) - 1)/1E4
  A <- sweep(A, 1, og_row_sums, '*')
  return(A)
}