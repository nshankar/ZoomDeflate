# Nikhil's wd
setwd('/Users/nikhil/Documents/College/Math 651/ZoomDeflate/')

library(reactable)

# Data sets to run code upon
nGroups <- c(2, 5, 10)
nCells <- c(1000, 10000)
nGenes <- c(5000, 1000)

data <- matrix(NA, nrow = 6, ncol = 5)
row_names <- rep("", 6)
ARI_col_order <- c(1,3,5,4,2)

for (j in 1:length(nCells)) {
  for (i in 1:length(nGroups)) {
    ARI <- read.csv(paste("Clustering_Data_dim_100/", nGroups[i], "_groups_", 
                              nCells[j], "_cells_", nGenes[j], "_genes/ARI_report_csv", sep=""))
    ID = paste("(", nGroups[i], ", ", nCells[j], ", ", nGenes[j], ")", sep="")
    row = i + (j-1)*length(nGroups)
    row_names[row] <- ID
    data[row, ] <- ARI[,5][ARI_col_order]
  }
}

data <- as.data.frame(round(data, 2), row.names=row_names)

orange_pal <- function(x) rgb(colorRamp(c("ivory", "#c9ecb4"))(x), maxColorValue = 255)
maxWidth_ = 90
mystyle <- function(value) {
  normalized <- ((value - min(data)) / (max(data) - min(data)))
  color <- orange_pal(normalized)
  list(background = color)
}

reactable(
  data, rownames=TRUE,
  style = list(fontFamily = "Work Sans, sans-serif", fontSize = "14px"),
  columns = list(
    .rownames = colDef(name = "(# Groups,    # Cells, # Genes)",
                       headerStyle = list(fontSize = "12px")),
    V1 = colDef(name = "Truth", style = mystyle, 
                maxWidth = maxWidth_, headerStyle = list(fontSize = "12px")),
    V2 = colDef(name = "ZoomDeflate", style = mystyle, 
                maxWidth = maxWidth_, headerStyle = list(fontSize = "12px")),
    V3 = colDef(name = "ALRA-Alt", style = mystyle, 
                maxWidth = maxWidth_, headerStyle = list(fontSize = "12px")),
    V4 = colDef(name = "ALRA", style = mystyle, 
                maxWidth = maxWidth_, headerStyle = list(fontSize = "12px")),
    V5 = colDef(name = "Observed", style = mystyle, 
                maxWidth = maxWidth_, headerStyle = list(fontSize = "12px"))
  ),
  columnGroups = list(
    colGroup(name = "Adjusted Rand Indices (10 k-means)", columns = c("V1", "V2", "V3", "V4", "V5"), align = "center")
  )
)
