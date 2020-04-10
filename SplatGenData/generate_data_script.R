library(splatter)

### Parameters (choose wisely) ###
random_seed <- 43
num_cells <- 10000
num_genes <- 1000
num_groups <- 10
dropout <- 1 #boolean valued for now
dropout.mid.value <- 2
dropout.shape.value <- -1

PATHNAME <- "/Users/nikhil/Documents/College/Math 651/ZoomDeflate/SplatGenData/10_groups_10000_cells_1000_genes/"

### Run simulation ###
params <- newSplatParams()
params <- setParam(params, "seed", random_seed)
params <- setParam(params, "batchCells", num_cells)
params <- setParam(params, "nGenes", num_genes)
if(num_groups > 1) {
  # Uniformly distributed groups
  probabilities = rep(1/num_groups, num_groups)
  params <- setParam(params, "group.prob", probabilities)
  if (dropout) {
    dropout.mid = rep(dropout.mid.value, num_groups) #default, possibly a bad assumption
    dropout.shape = rep(dropout.shape.value, num_groups) #default, possibly a bad assumption
    params <- setParam(params, "dropout.mid", dropout.mid)
    params <- setParam(params, "dropout.shape", dropout.shape)
    params <- setParam(params, "dropout.type", "group")
  }
  sim <- splatSimulate(params, method = "group")
} else {
  if (dropout) {
    params <- setParam(params, "dropout.mid", dropout.mid.value)
    params <- setParam(params, "dropout.shape", dropout.shape.value)
    params <-setParam(params, "dropout.type", "experiment")
  }
  sim <- splatSimulate(params)
}

sum(rowSums(counts(sim) == 0))

### Save data and metadata (which metadata is important?) ###
dir.create(PATHNAME)
write.csv(counts(sim), paste(PATHNAME, "counts.csv"))
if(dropout) {
  write.csv(assays(sim)[["TrueCounts"]], paste(PATHNAME, "true_counts.csv"))
  write.csv(assays(sim)[["Dropout"]]*1, paste(PATHNAME, "dropouts.csv"))
}

