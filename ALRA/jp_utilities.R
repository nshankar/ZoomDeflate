# Inputs: mask (logical matrix with TRUE in locations of dropouts)
#         truth (non-negative valued matrix with true counts. Should have same scaling as recon)
#         recon (non-netative valued matrix with reconstructed counts.)

# Outputs: f
zero_quality_stats <- function(mask, truth, recon) {
  # Calculate biological and technical zeros
  bio_zeros_mask <- (truth == 0)
  tech_zeros_mask <- (truth != 0) & (mask)
  
  frac_bio_zeros <- sum(bio_zeros_mask) / length(truth)
  frac_tech_zeros <- sum(tech_zeros_mask) / length(truth)

  # Count biological zeros preserved
  bio_zeros_preserved <- sum(A_norm_completed[bio_zeros_mask] == 0) 
  frac_bio_zeros_preserved <- bio_zeros_preserved / sum(bio_zeros_mask)

  # Count technical zeros reconstructed
  tech_zeros_reconned <- sum(A_norm_completed[tech_zeros_mask] != 0) 
  frac_tech_zeros_reconned <- tech_zeros_reconned / sum(tech_zeros_mask)

return (list(frac_bio_zeros_preserved = frac_bio_zeros_preserved,
             frac_tech_zeros_reconned = frac_tech_zeros_reconned))
 }
 
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