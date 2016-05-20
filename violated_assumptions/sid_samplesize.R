library(bnlearn)
library(pcalg)
library(igraph)
library(graph)
library(SID)

#source("test_utils.R")

bn.get.sid = function(fit, true.am, method.name){
  #fit: a fitted object from bnlearn
  #true.am: adjacency matrix of ground truth 
  #method.name: str name of method used in fitting the object fit
  #return: sid object from structIntervDist estimate to true.am
  
  cpdag = method.name == "hc" # the hill climbing algorithm returns a DAG
                              # estimation without interventoinal data can not be better than a CPDAG
                              # thus the estimate is converted to the corresponding CPDAG and then analysed
  
  am = get.adjm.from.bn(fit, true.am, cpdag) #@test_utils.R
  sid= structIntervDist(trueGraph = true.am, estGraph = am)
  return(sid)
}

pcalg.get.sid = function(fit, true.am, method.name){
  #fit: a fitted object from pcalg
  #true.am: adjacency matrix of ground truth 
  #method.name: str name of method used in fitting the object fit
  #return: sid object from structIntervDist estimate to true.am
  if(method.name == "pc"){
    am = get.adjacency(igraph.from.graphNEL(fit@graph))
  }else if(method.name == "fci"){
    am = fit@amat
  }else if( method.name %in% c("ges", "gies")){
    am = get.adjm.from.essgraph(fit$essgraph)
  } 
  sid = structIntervDist(trueGraph = true.am, estGraph = am)
  
  return(sid)

}

sid.size = function(data, method, method.name, true.am){
  #data: data frame with data for estimation
  #method: variable of method in bnlearn or pcalg
  #method.name: str name of method
  #true.am: adjacency matrix of ground truth
  if(method.name %in% c("gs", "iamb", "fast.iamb", "inter.iamb", "hc")){
    fit = method(data)
    sid = bn.get.sid(fit, true.am, method.name)
  }
  else if(method.name %in% c("pc", "ges", "fci")){
    fit = fit.pcalg(data, method, method.name)
    sid = pcalg.get.sid(fit, true.am, method.name)
  }
  
  return(sid)
}

get.upper = function(sid){
  return(sid$sidUpperBound)
}
get.lower = function(sid){
  return(sid$sidLowerBound)
}
get.sid = function(sid){
  return(sid$sid)
}

get.eq.upper = function(true.am){
  #true.am: adjacency matrix of a DAG
  #return: max sid of a dag in the true eq class to ground truth
  g = graphAM(true.am, edgemode = "directed")
  g = as(g, "graphNEL")
  am = dag2cpdag(g)
  
  return(structIntervDist(true.am, am)$sidUpperBound)
}

plot.sid= function(dat.list, size, method, method.name, true.am){
  #dat.list list of data frames of samples of different size from same sem
  #size: vector of nrow of data frames in dat.list
  #method: variable of method in bnlearn or pcalg
  #method.name: str 
  #true.am: adjacency matrix of sem
  
  sid = lapply(dat.list, sid.size, method = method, method.name = method.name, true.am = true.am)
  upper = sapply(sid, get.upper)
  lower = sapply(sid, get.lower)
  #es.sid = sapply(sid, get.sid)
  eq.upper.lim = get.eq.upper(true.am)
  ymax = max(max(upper), eq.upper.lim)
  
  plot(size, upper, pch = 20, cex = 0.4, main = method.name, 
       ylim = c(0, ymax), xlab = "sample size", col = "red",
       ylab = "SID lower, upper bound")
  points(size, lower, pch = 20, cex = 0.4, col = "green4")
  abline(h=eq.upper.lim, lty = 2)
}

sid.int = function(sample, data, target.index, int.targets, method, true.am){

  score = new("GaussL0penIntScore", as.matrix(data[[sample]]), targets = int.targets, target.index = target.index[[sample]])
  fit = method(ncol(data[[sample]]), int.targets, score)
  am = get.adjm.from.essgraph(fit$essgraph)
  
  sid= structIntervDist(trueGraph = true.am, estGraph = am)
  
  return(sid)
}

plot.int.sid = function(dat.list, int.index.list, int.targets, method, method.name, true.am, size){
  
  
  sample = 1:length(dat.list)
  sid = lapply(sample, FUN = sid.int, data = dat.list, target.index = int.index.list, int.targets = int.targets, method = method,  true.am = true.am)

  upper = sapply(sid, get.upper)
  lower = sapply(sid, get.lower)
  #es.sid = sapply(sid, get.sid)
  num.interventions = length(int.targets)
  
  eq.upper.lim = get.eq.upper(true.am)
  ymax = max(max(upper), eq.upper.lim)
  
  plot(num.interventions*size, upper, pch = 20, cex = 0.4, main = method.name, 
       ylim = c(0, ymax), xlab = "sample size", col = "red",
       ylab = "SID lower, upper bound")
  points(num.interventions*size, lower, pch = 20, cex = 0.4, col = "green4")
  abline(h=eq.upper.lim, lty = 2)
  
}


