# Nikhil's wd
# setwd('/Users/nikhil/Documents/College/Math 651/ZoomDeflate/softImputeRecon')

# Jeremy's wd
setwd('~/Documents/Projects/ZoomDeflate/')

library(softImpute)
library(ggplot2)
library(dplyr)
library(hrbrthemes)
source('ALRA/alra.R')
source('ALRA/jp_utilities.R')
set.seed(41) # Some ML projects make this seed a hyper-parameter

data_path = "/SplatGenData/5_groups_1000_cells_5000_genes/"
MASK_PATHNAME = paste(data_path,"dropouts.csv", sep="")
DATA_PATHNAME = paste(data_path,"counts.csv", sep="")
TRUE_PATHNAME = paste(data_path,"true_counts.csv",sep="")

# load matrices in format: rows are genes, columns are cells
mask <- as.matrix(read.csv(MASK_PATHNAME, header=FALSE, sep=" "))
data <- as.matrix(read.csv(DATA_PATHNAME, header=FALSE, sep=" "))
truth <- as.matrix(read.csv(TRUE_PATHNAME, header=FALSE, sep=" "))

### RUN ALRA ### 


### RUN SOFTIMPUTE ### 

