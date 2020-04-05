library(splatter)

### Parameters (choose wisely) ###
random_seed <- 43
num_cells <- 1000
num_genes <- 50
num_groups <- 2
dropout <- 1 #boolean valued for now

PATHNAME <- "/Users/nikhil/Documents/College/Math 651/generate_data/example_data/"

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
    dropout.mid = rep(0, num_groups) #default, possibly a bad assumption
    dropout.shape = rep(-1, num_groups) #default, possibly a bad assumption
    params <- setParam(params, "dropout.mid", dropout.mid)
    params <- setParam(params, "dropout.shape", dropout.shape)
    params <- setParam(params, "dropout.type", "group")
  }
  sim <- splatSimulate(params, method = "group")
} else {
  if (dropout) {
    params <-setParam(params, "dropout.type", "experiment")
  }
  sim <- splatSimulate(params)
}

### Save data and metadata (which metadata is important?) ###
dir.create(PATHNAME)
write.csv(counts(sim), paste(PATHNAME, "counts.csv"))
if(dropout) {
  write.csv(assays(sim)[["TrueCounts"]], paste(PATHNAME, "true_counts.csv"))
  write.csv(assays(sim)[["Dropout"]], paste(PATHNAME, "dropouts.csv"))
}

