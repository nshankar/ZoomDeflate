library(splatter)

### Parameters (choose wisely) ###
random_seed <- 43
num_cells <- 1000
num_genes <- 10000
num_groups <- 10
dropout <- 1 #boolean valued for now
dropout.mid.value <- 2
dropout.shape.value <- -1

PATHNAME <- paste("/Users/nikhil/Documents/College/Math 651/ZoomDeflate/SplatGenData/", 
                  num_groups,"_groups_", num_cells, "_cells_", num_genes, "_genes/", sep="")
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
write.table(counts(sim), paste(PATHNAME, "counts.csv", sep=""), row.names=FALSE, col.names=FALSE)
if(dropout) {
  write.table(assays(sim)[["TrueCounts"]], paste(PATHNAME, "true_counts.csv", sep=""), row.names=FALSE, col.names=FALSE)
  write.table(assays(sim)[["Dropout"]]*1, paste(PATHNAME, "dropouts.csv", sep=""), row.names=FALSE, col.names=FALSE)
}
if(num_groups > 1) {
  groupData = sapply(sim@colData@listData[["Group"]], function(s) 
                  for (i in 1:num_groups) {
                    if (s == paste("Group", i, sep="")) {
                      return(i)
                    }
                  }
                )
                  
  write.table(groupData, paste(PATHNAME, "group_data.csv", sep=""), row.names=FALSE, col.names=FALSE)
}

