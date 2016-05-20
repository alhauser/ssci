library(igraph)
source("../../test_utils.R")
mmat = as.matrix(read.table("mmat.txt"))

adjm = adjm.from.mmat(mmat)

logistic = t(as.matrix(read.table("logistic.txt")))

g = graph.adjacency(adjm, mode = "directed")

col = rep("black", ecount(g))
edge.nr = 1
for(i in 1:nrow(logistic)){
  for(j in 1:ncol(logistic)){
    if(logistic[i,j] == 1){
      col[edge.nr] = "red"
      edge.nr = edge.nr + 1
    }else if(adjm[i,j] == 1){
      edge.nr = edge.nr + 1
    }
  }
}

pdf("final_graph.pdf")
plot.igraph(g, edge.color = col)
dev.off()
