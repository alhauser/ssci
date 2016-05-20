#####################################
# functions to estimate different scores of an inferred model
# requires some functions of test_utils.R
#########################################


#source("../test_utils.R")


library(CompareCausalNetworks)
library(bnlearn)
library(igraph)
library(graph)
library(pcalg)
library(SID)

# -----------------------------------------------------
# functions for sub- and resampling
# -----------------------------------------------------
get.sampl = function(sampl.index, data){
  # sample.index: vector of indices of subsample
  # data: data frame of full sample
  return(data[sampl.index,])
}

subsample = function(n, dil.fact, data){
  # data: the data set to sub sample from
  # dil.fact : dilution factro, prob for a sample to be in the sub sample
  # n : number of subsamples
  # return: list of n subsamples
  
  sub.samp.size = round(nrow(data)*dil.fact)
  sub.index = replicate(n, sample(1:nrow(data), sub.samp.size, replace = F), simplify = F)# generate n subsets of indeices without replacement
  subsamp = lapply(sub.index, get.sampl, data)
  return(subsamp)
}

resample = function(n, data, resamp.size){
  # n: number of resampless
  # resamp.size : number of samples in the resample
  # return: a list of n resampled data sets
  size = nrow(data)
  resamp.index = replicate(n, sample(1:size, resamp.size, replace = T), simplify = F)# generate n resamples of the indices with replacement
  #print(resamp.index)
  resamp = lapply(resamp.index, get.sampl, data)
  return(resamp)
}

# -----------------------------------------------------
# functions for sub- and resamping of interventional data
# -----------------------------------------------------

concatenate.experiments = function(experiment.list){
  
  #experiment.list: list with experiments (list( df, interventions))
  #return: list of 2, data frame with data from all experiments in experiment.list, 
  # list of nrow(data frame) lists with vector of intrerventions of every data point in data frame
  
  #   print(experiment.list[[1]])
  #   str(experiment.list)
  df = data.frame(replicate(ncol(experiment.list[[1]][[1]]), numeric(0)))
  
  intervention.list = list()
  for(experiment in experiment.list){
    df = rbind(df, experiment[[1]])
    intervention.list = append(intervention.list, replicate(nrow(experiment[[1]]), experiment[[2]], simplify = F))    
  }
  return(list(df, intervention.list))
}

concatenate.experiments2 = function(experiment.list){
  #experiment.list: list with experiments (list( df, interventions))
  #return: list of 2, data frame with data from all experiments in experiment.list, 
  # named lists with vectors of indices of datapoints in df with intervention on variable with corresponding name
  # formated for tabu search in bnlearn
  
  df = data.frame(replicate(ncol(experiment.list[[1]][[1]]), numeric(0)))
  
  intervention.list = list()
  int.names = numeric(0)
  pos = 0
  for(experiment in experiment.list){
    df = rbind(df, experiment[[1]])
    if(experiment[[2]] != 0){ # experiment has intervention
      # append list vector of position of data in complete data set to intervention.list
      intervention.list = append(intervention.list, list((pos+1):(pos+nrow(experiment[[1]]))))
      # get name of targetnode
      int.names = append(int.names, names(experiment[[1]])[experiment[[2]]])
    }
  pos = nrow(df)
  }
  
  names(intervention.list) = int.names
  
  return(list(df, intervention.list))
  
  
}

subsample.experiment.list = function(experiment.list, dil.fact = 0.5){
  
  #experiment.list: list with experiments (list( df, interventions))
  #dil.fact: float [0,1] fraction of samplesize of original data set
  #return: list  of two 1: data frame with subsample, 2: list of interventions
  subsample.list = list()
  c = 1
  for(experiment in experiment.list){
    #str(experiment)
    subsample.list[[c]] = list(subsample(1, dil.fact, experiment[[1]])[[1]], experiment[[2]])
    c = c + 1
  }
  
  return(concatenate.experiments(subsample.list))
  
}

resample.experiment.list = function(experiment.list){
  #experiment.list: list with experiments (list( df, interventions))
  #return: list two 1: data frame with subsample, 2: list of interventions
  
  resamp.list = list()
  c = 1
  for(experiment in experiment.list){
    resamp.list[[c]] = list(resample(1, experiment[[1]], nrow(experiment[[1]]))[[1]], experiment[[2]])
    c = c + 1
  }
  
  return(concatenate.experiments(resamp.list)) 
}

#######################################################################
# -----------------------------------------------------
# scoring functions
# -----------------------------------------------------

# -----------------------------------------------------
# scoring functions based on shd
# -----------------------------------------------------
struct.ham.dist = function(combination, adj.list){
  # combinations: vector of two index in adj.list of adjacency matixes to compare
  # adj.list: list of adjacency matices
  am1 = adj.list[[combination[1]]]
  am2 = adj.list[[combination[2]]]
  
  return(sum(am1 != am2))
}

pairwise.struct.ham.dist = function(adj.list){
  #adj.list: list of estimated adjacency matrices 
  #return: vector of pairvise distances
  combinations = combn(length(adj.list), 2)
  struct.dist = apply(combinations, 2 , struct.ham.dist, adj.list)
}


struct.ham.dist.stability = function(dat, resamp.method, estimation.method,
                     n = 100, dil.fact = 0.9, resamp.size = nrow(dat)){
  # dat: data frame of the data to estimate the score
  # resamp.method: str subsampling(without replacement) or resampling(with replacement)
  # estimation.mathod: str implementation of an algorirthm in bnlearn or CompareCausalNetworks
  # n: int number of subsamples, resamples respectively
  # dil.fact: numeric [0,1] dilution for subsampling
  # resamp.size: size of resampled data sets
  # return : average pairwise distance per possible edge
  samp.list = switch(resamp.method,
                     resample = {resample(n, dat, resamp.size)},
                     subsample = {subsample(n, dil.fact, dat)},
                     stop("No valid resample method given!"))
  
  adj.list = lapply(samp.list, est.adjm, estimation.method) #est.adjm: function @ test_utils
  dist = pairwise.struct.ham.dist(adj.list)
  return(mean(dist)/(ncol(dat)^2))
}

# -----------------------------------------------------
# scoring functions base on StARS
# -----------------------------------------------------

stars.score = function(adj.list){
  #adj.list: list of adjacency matrix
  #return: score simillar to the one used in StARS
  
  # find probability of edge to be selected
  prob.mat = matrix(0, nrow = nrow(adj.list[[1]]), ncol = ncol(adj.list[[1]]))
  for(i in adj.list){
    prob.mat = prob.mat + i
  }
  prob.mat = prob.mat/length(adj.list)
  
  score.mat = 2*prob.mat*(1-prob.mat) # probability for an edge to be different in an other estimation
  sscore = sum(score.mat)/(ncol(prob.mat)^2 - ncol(prob.mat))# averaging over all possible edges in directed graph
                                                  # edge to one self not possible thus ncol(p.m)^2 - ncol(p.m) possible edges
  
  return(sscore)
}

est.stars.score = function(dat, resamp.method = "subsample", estimation.method,
                           n = 200, dil.fact = 0.5, resamp.size = nrow(dat)){
  # dat: data frame of the data to estimate the score
  # resamp.method: str subsampling(without replacement) or resampling(with replacement)
  # estimation.mathod: str implementation of an algorirthm in bnlearn or CompareCausalNetworks
  # n: int number of subsamples, resamples respectively
  # dil.fact: numeric [0,1] dilution for subsampling
  # resamp.size: size of resampled data sets
  # return : score simillar to the one used in StARS
  
  samp.list = switch(resamp.method,
                     resample = {resample(n, dat, resamp.size)},
                     subsample = {subsample(n, dil.fact, dat)},
                     stop("No valid resample method given!"))
  adj.list = lapply(samp.list, est.adjm, estimation.method) #est.adjm: function @ test_utils
  score = stars.score(adj.list)
  return(score)
  
}


num.edges = function(adjm.list){
  # adj.list: list of adjacency matrixes
  # return: vector of number of edges 
  return(sapply(adjm.list, function(x) sum(x)))
}

score.edgedist = function(dat, resamp.method, estimation.method,
                          n = 200, dil.fact = 0.5, resamp.size = nrow(dat)){
  
  # dat: data frame data used in analysis
  # resamp.method: str "subsample" "resample"
  # n: int number of sub/re samples
  # dil.fact: float [0,1] dilution in subsampling
  # resamp.size: int size of resamples
  # return : list of 2 stars score, average directed density of estimates used for score
  
  samp.list = switch(resamp.method,
                     resample = {resample(n, dat, resamp.size)},
                     subsample = {subsample(n, dil.fact, dat)},
                     stop("No valid resample method given!"))
  adj.list = lapply(samp.list, est.adjm, estimation.method) #est.adjm: function @ test_utils
  score = stars.score(adj.list)
  density = mean(sapply(adj.list, directed.density))
  return(list(score, density))
}


directed.density = function(adjm){
  # adjm: adjacency matrix
  # return: directed density ( number of edges per possible edges, undirected edges are only counted as one)
  edges = sum(adjm)
  undir.edges = 0
  for(i in 1:nrow(adjm)){
    for(j in 1:ncol(adjm)){
      if(i < j && adjm[i,j] == 1 && adjm[j,i] == 1){
        undir.edges = undir.edges + 1
      }
    }
  }
  
  #print(undir.edges)
  return((edges - undir.edges ))#/choose(ncol(adjm), 2))
}

scaled.score = function(dat, resamp.method, estimation.method,
                        n = 200, dil.fact = 0.5, resamp.size = nrow(dat)){
  
  # dat: data frame data used in analysis
  # resamp.method: str "subsample" "resample"
  # n: int number of sub/re samples
  # dil.fact: float [0,1] dilution in subsampling
  # resamp.size: int size of resamples
  # return : scaled score stars/stars of random dags of same density
  samp.list = switch(resamp.method,
                     resample = {resample(n, dat, resamp.size)},
                     subsample = {subsample(n, dil.fact, dat)},
                     stop("No valid resample method given!"))
  adj.list = lapply(samp.list, est.adjm, estimation.method) #est.adjm: function @ test_utils
  
  score = stars.score(adj.list)
  nnodes = ncol(dat)
  est.density = mean(sapply(adj.list, directed.density))
  #print(nnodes)
  #print(est.density)
  random.adj.list = replicate(n, randomDAG(nnodes, est.density), simplify = F)
  random.score = stars.score(random.adj.list) 
  return(score / random.score )
}



# -----------------------------------------------------
# functions for stars scores on interventional data
# -----------------------------------------------------
stars.intervention.list = function(index, data.list, intervention.list, 
                                   resamp.method, n = 200, dil.fact = 0.5, resamp.size = nrow(dat)){
  return(stars.intervention.edgelist(data.list[[index]], intervention.list[[index]], resamp.method, n, dil.fact, resamp.size))
}


fit.from.data.list = function(samp.list, method){
  # str(samp.list)
  getParents(samp.list[[1]], method = method, interventions = samp.list[[2]])
}
stars.intervention.edgelist = function(samp.list, method){

  adj.list = lapply(samp.list, fit.from.data.list, method)
  score = stars.score(adj.list)
  density = mean(sapply(adj.list, directed.density))
  return(list(score, density))
  
}

# -----------------------------------------------------
# scoring functions based on sid
# -----------------------------------------------------
sym.sid= function(combination, adj.list){
  
  # combinations: vector of two index in adj.list of adjacency matixes to compare
  # adj.list: list of adjacency matices
  # return: symetric average of uppwer and lower bounds of sid
  
  adjm1 = adj.list[[combination[1]]]
  adjm2 = adj.list[[combination[2]]]
  sid1 = structIntervDist(adjm1, adjm2)
  sid2 = structIntervDist(adjm2, adjm1)
  return((sid1$sidUpperBound + sid1$sidLowerBound + sid2$sidUpperBound +sid2$sidLowerBound)/4)
}
  

pairwise.sym.sid = function(adj.list){
  # adj.list: list of adjacency matrixes
  # return: average symetized structural intervention distance of all pairwise comparisons
  
  combinations = combn(length(adj.list), 2)
  struct.dist = apply(combinations, 2, sym.sid, adj.list)
  return(mean(struct.dist))
}

est.pw.sym.sid = function(dat, resamp.method, estimation.method,
                          n = 100, dil.fact = 0.5, resamp.size = nrow(dat)){
  # dat: data frame of the data to estimate the score
  # resamp.method: str subsampling(without replacement) or resampling(with replacement)
  # estimation.mathod: str implementation of an algorirthm in bnlearn or CompareCausalNetworks
  # n: int number of subsamples, resamples respectively
  # dil.fact: numeric [0,1] dilution for subsampling
  # resamp.size: size of resampled data sets
  # return: average symetized structural intervention distance of all pairwise comparisons of 
  #         estimates on subsamples
  
  samp.list = switch(resamp.method,
                     resample = {resample(n, dat, resamp.size)},
                     subsample = {subsample(n, dil.fact, dat)},
                     stop("No valid resample method given!"))
  adj.list = lapply(samp.list, est.adjm, estimation.method) #est.adjm: function @ test_utils
  pw.dist = pairwise.sym.sid(adj.list)
  return(pw.dist)
  
}
