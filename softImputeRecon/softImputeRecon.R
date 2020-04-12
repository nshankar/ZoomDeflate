require(softImpute)

# Currently cheating by using the true dropouts as a mask
# In the future, use predicted dropouts as a mask
### PATHNAME is currently for Nikhil's computer, fix soon ###
PATHNAME = "/Users/nikhil/Documents/College/Math 651/ZoomDeflate/SplatGenData/5_groups_1000_cells_5000_cells/"
MASK_PATHNAME = paste(PATHNAME, "dropouts.csv", sep="")
DATA_PATHNAME = paste(PATHNAME, "counts.csv", sep="")


set.seed(1011)# Some ML projects make this seed a hyper-parameter
lambda_= 1.9  # regularization weight emphasizing the nuclear norm
type_ = "als" # either als or svd, als is faster
rank.max_ = 0 # not a necessary argument, 
              # but an estimate could be useful

normalize_rows_by_noise = 0 #boolean


mask <- as.matrix(read.csv(MASK_PATHNAME, row.names=1))
data <- as.matrix(read.csv(DATA_PATHNAME, row.names=1))

if(normalize_row_by_noise) {
  # todo
}

og_data <- data # for later comparison
data[as.logical(mask)] = NA

if (rank.max_) {
  fits = softImpute(data, rank.max = rank.max_,
                    lambda = lambda_, trace=TRUE, type=type_)
} else {
  fits = softImpute(data, lambda = lambda_,
                    trace=TRUE, type=type_)
}

output <- complete(data, fits)