library(plyr)
library(graph)
library(pcalg)
library(SID)
source("../../test_utils.R")

set.seed(6)
n= 500
rand.dags = rlply(n, as.matrix(randomDAG(10, 0.3)))
rand.cpdags = rlply(n, cpdag.from.adjm(as.matrix(randomDAG(10,0.3))))

densities = runif(n, 0.05, 0.95)
rand.dags = lapply(densites, randomDAG, p= 10)
rand.dags = lapply(rand.dags, as.matrix)
rand.cpdags =  lapply(densites, randomDAG, p= 10)
rand.cpdags = lapply(rand.cpdags, as.matrix)
rand.cpdags = lapply(rand.cpdags, cpdag.from.adjm)
# pw.eq.shd = function(adjm.list){
#  
#   eq.shd = numeric(length(adjm.list)^2)
#   c = 0
#   for(i in adjm.list){
#    for (j in adjm.list){
#      cpdag = cpdag.from.adjm(i)
#      cpdag = as(cpdag, "graphAM")@adjMat
#      eq.shd[c] = sum( cpdag != j)
#      c = c+1
#    }
#   } 
#   return(eq.shd)
# }
# 
# pw.eq.distance = function(adjm.list){
#   
#   eq.dist = numeric(length(adjm.list)^2)
#   c = 0
#   for(i in adjm.list){
#     for (j in adjm.list){
#       cpdag = cpdag.from.adjm(i)
#       cpdag = as(cpdag, "graphAM")@adjMat
#       sid = structIntervDist(cpdag, j)
#       upper =
#       eq.dist[c] = 
#       c = c+1
#     }
#   } 
#   return(eq.dist)
#   
# }

est.eqclass.sid = function(index, dag, cpdag){
  sid = structIntervDist(dag[[index]], cpdag[[index]])
  true.upper = true.eq.upper(dag[[index]])
  dist = sid$sidLowerBound + abs(sid$sidUpperBound - true.upper)
  diff = sid$sidLowerBound + sid$sidUpperBound - true.upper
  
  return(list(dist, diff))
}

est.eqclass.diff = function(index, dag, cpdag){
  sid = structIntervDist(dag[[index]], cpdag[[index]])
  true.upper = true.eq.upper(dag[[index]])
  diff = sid$sidLowerBound + sid$sidUpperBound - true.upper
  return(diff)
}


est.shd = function(index, dag, cpdag){
  shd = sum(dag[[index]] != cpdag[[index]])
  return(shd)
}
eqclass.shd = function(index, dag, cpdag){
  true.cpdag = cpdag.from.adjm(dag[[index]])
  shd = sum(true.cpdag != cpdag[[index]])
  return(shd)
}

sid = sapply(1:n, est.eqclass.sid, rand.dags, rand.cpdags)
sid.diff = sid[[2]]
sid.dist = sid[[1]]

shd = sapply(1:n, est.shd, rand.dags, rand.cpdags)
shd.eqclass= sapply(1:n, eqclass.shd, rand.dags, rand.cpdags)

scores = data.frame(shd, shd.eqclass, sid.diff, sid.dist )
pdf("scores.pdf")
pairs(scores, pch = 20, cex = 0.7)
dev.off()
