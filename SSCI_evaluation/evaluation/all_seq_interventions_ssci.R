source("../../test_utils.R")
source ("../../scoring.R")

library(SID)

est.shd = function(index, est.adjm.list, true.adjm.list){
  shd = sum(est.adjm.list[[index]] != true.adjm.list[[index]])
  return(shd)
}
eqclass.shd = function(index, est.adjm.list, true.adjm.list){
  cpdag = cpdag.from.adjm(true.adjm.list[[index]])
  #cpdag = as(cpdag, "graphAM")@adjMat
  #print(cpdag)
  shd = sum(est.adjm.list[[index]] != cpdag)
  return(shd)
}
# eqclass shd is a bad measure like this, for sure for gies, because the true dag has a larger value than the 
# eqclass cpdag it self gies can get better than just the cpdag!

eq.to.eqclass.shd = function(index, est.adjm.list, true.adjm.list){
  true.cpdag = cpdag.from.adjm(true.adjm.list[[index]])
  #cpdag = as(cpdag, "graphAM")@adjMat
  est.cpdag = cpdag.from.adjm(as.matrix(est.adjm.list[[index]]))
  shd = sum(est.cpdag != true.cpdag)
  return(shd)
}
est.sid.measures = function(index, est.adjm.list, true.adjm.list){
  sid = structIntervDist(true.adjm.list[[index]], est.adjm.list[[index]])
  true.upper = true.eq.upper(true.adjm.list[[index]])
  dist = sid$sidLowerBound + abs(sid$sidUpperBound - true.upper)
  diff = sid$sidLowerBound + sid$sidUpperBound - true.upper
  return(list(dist, diff))
}

get.int.targets = function(n, nintervent){
  
  int.targets = sample(c(rep(1, nintervent), rep(0, n - nintervent)))
  
}

int.gaus.index = function(index, sem, int.targets, target.value = 1, effective.size.seq){
  
  node = 1
  d = data.frame(replicate(ncol(sem[[index]][[1]]), numeric(0)))
  for(i in int.targets[[index]]){
    if(i == 1){
      new.dat = get.int.lin.gaus(sem[[index]][[1]], node, target.value, sem[[index]][[2]], size = effective.size.seq[[index]])
      #str(new.dat)
     
    } else{
      new.dat = get.lin.gaus(sem = sem[[index]], size = effective.size.seq[[index]])
    }
    d = rbind(d, new.dat)
    node = node + 1
  }
  return(d)
}

int.exp.index = function(index, sem, int.targets, target.value = 1, effective.size.seq){
  
  node = 1
  d = data.frame(replicate(ncol(sem[[index]][[1]]), numeric(0)))
  for(i in int.targets[[index]]){
    if(i == 1){
      new.dat = get.int.exp(sem[[index]][[1]], node, target.value, sem[[index]][[2]], size = effective.size.seq[[index]])
      #str(new.dat)
      
    } else{
      new.dat = get.lin.exp(mmat = sem[[index]][[1]], lambda = 1 / sem[[index]][[2]], size = effective.size.seq[[index]])
    }
    d = rbind(d, new.dat)
    node = node + 1
  }
  
  return(d)
}
# c = 1
# node = 1
# for( i in int.targets[[1]]){
#   if(i == 1){
#     target.index[c:(c+4)] = rep(node, 5)
#   } else{
#     target.index[c:c+4] = rep(0, 5)
#   }
#   c = c +5
#   node = node+1
# }

get.intervention.list = function(index, int.targets, effective.size.seq){
  
  nnodes = length(int.targets[[index]])
  size = effective.size.seq[[index]]
  int.index = as.list(rep(int.targets[[index]] * (1:nnodes), rep(size,nnodes)))
  return(int.index)
}

est.adjm.gies = function(index, dat, intervention.list){
  return(getParents(dat[[index]], interventions = intervention.list[[index]], method = "gies"))
}

get.gaus.experiment = function(sem, interventions, tar.values, size){
  
  # sem: list of 2, matrix describing interactions in model, vector of standart deviations of nodes
  # interventions: vector of nodes to intervene on
  # target.values: vector with values for intervention nodes
  # size: int, sample size number of points to draw from distribution
  # return: list of 2, data frame of with sample, vector of interventions
  
  df = get.int.lin.gaus(sem[[1]], interventions, tar.values, sem[[2]], size)
  return(list(df, interventions))
}

get.exp.experiment = function(sem, interventions, tar.values, size){
  # sem: list of 2, matrix describing interactions in model, vector of standart deviations of nodes
  # interventions: vector of nodes to intervene on
  # target.values: vector with values for intervention nodes
  # size: int, sample size number of points to draw from distribution
  # return: list of 2, data frame of with sample, vector of interventions
  
  df = get.int.exp(sem[[1]], interventions, tar.values, 1/sem[[2]], size)
  return(list(df, interventions))
}

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


subsample.experiment.list = function(experiment.list, dil.fact){
  
  #experiment.list: list with experiments (list( df, interventions))
  #n : int, number of subsamples
  #dil.fact: float [0,1] fraction of samplesize of original data set
  #return: list of n lists with 1: subsample, 2: vectors of interventions
  subsample.list = list()
  c = 1
  for(experiment in experiment.list){
    subsample.list[[c]] = list(subsample(1, dil.fact, experiment[[1]])[[1]], experiment[[2]])
    c = c + 1
  }
  
  return(concatenate.experiments(subsample.list))
  
}

resample.experiment.list = function(experiment.list){
  
  resamp.list = list()
  c = 1
  for(experiment in experiment.list){
    resamp.list[[c]] = list(resample(1, experiment[[1]], nrow(experiment[[1]]))[[1]], experiment[[2]])
    c = c + 1
  }
  
  return(concatenate.experiments(resamp.list)) 
}
score.intervention =  function(nnodes, nintervent, target.value,
                               densities, size.seq, effect.rnge, noise.rnge, method, resamps, dil.fact){
  #n: int number of variables in sem, nodes in graph
  #density: R[0,1] density of graph, fraction of edges to number of possible edges
  #effect.rnge: vector of upper and lower limit to sample magnitude of causal effects from
  #noise.rnge: vector of upper and lower limit oto sample the standardvariation of the noise terms
  #trails: int numer of different sem to test
  #method: str name of algo to infere model with ( CompareCausalNetworks or bnlearn)
  #resamp.method: str "subsample" or "resample"
  #resamps: int number of resamples used for score estimation
  #size: int number of samples drawn from every sem
  
  
  gaus.sem.list = all.seq.sem(nnodes, densities, size.seq, effect.rnge, noise.rnge)
  exp.sem.list = all.seq.sem(nnodes, densities, size.seq, effect.rnge, noise.rnge)
  mmat.list = lapply(append(gaus.sem.list, exp.sem.list), function(x) x[[1]])
  #sdev = lapply(sem, function(x) x[[2]])
  adjm.list = lapply(mmat.list, adjm.from.mmat)
  
  #   perm.list = replicate(trials, sample(1:n, replace = F), simplify = F)
  #   adjm.list = lapply(1:trials, scramble.adjmat.list, perm.list, adjm.list)
  #   
  
  effective.size.seq = rep(size.seq, length(nnodes)*length(densities))
  trials = length(effective.size.seq)
  
  gaus.int.targets = replicate(trials, c(0, sample(1:nnodes, nintervent)), simplify = F)
  
#   print(gaus.int.targets[[1]])
  exp.int.targets = replicate(trials, c(0, sample(1:nnodes, nintervent)), simplify = F)
  
  est.am.list = list()
  subs.res = list()
  resamp.res = list()
  for(i in 1:trials){
      gaus.experiment.list = list()
      
      gaus.subsamples.list = list()
      gaus.resamples.list = list()
      
      exp.experiment.list = list()
      
      exp.subsamples.list = list()
      exp.resamples.list = list()
      
      
      c = 1
      for(target in gaus.int.targets[[i]]){ 
        
        gaus.experiment.list[[c]] = get.gaus.experiment(gaus.sem.list[[i]], target, target.value, effective.size.seq[[i]])
        c = c + 1
      }
      c = 1
      for(target in exp.int.targets[[i]]){ 
        
        exp.experiment.list[[c]] = get.exp.experiment(exp.sem.list[[i]], target, target.value, effective.size.seq[[i]])
        c = c + 1
      }
      #print(experiment.list)
      gaus.full.data.set = concatenate.experiments(gaus.experiment.list)
#      str(gaus.full.data.set[[2]])
#        str(gaus.full.data.set[[1]])
      exp.full.data.set = concatenate.experiments(exp.experiment.list)
      
      est.am.list[[i]] = getParents(as.data.frame(gaus.full.data.set[[1]]), method = "gies", interventions = gaus.full.data.set[[2]], directed = F)
      est.am.list[[i + trials]] = getParents(as.data.frame(exp.full.data.set[[1]]),  method = "gies", interventions = exp.full.data.set[[2]], directed = F)
      
      gaus.subsamples.list = replicate(resamps, subsample.experiment.list(gaus.experiment.list, dil.fact), simplify = F)
      gaus.resample.list = replicate(resamps, resample.experiment.list(gaus.experiment.list), simplify = F)
      
      exp.subsamples.list = replicate(resamps, subsample.experiment.list(exp.experiment.list, dil.fact), simplify = F)
      exp.resample.list = replicate(resamps, resample.experiment.list(exp.experiment.list), simplify = F)
      
      subs.res[[i]] = stars.intervention.edgelist(gaus.subsamples.list, method = "gies")
      resamp.res[[i]] = stars.intervention.edgelist(gaus.resample.list, method = "gies")
      
      subs.res[[i + trials]] = stars.intervention.edgelist(exp.subsamples.list, method = "gies")
      resamp.res[[i + trials]] = stars.intervention.edgelist(exp.resample.list, method = "gies")
  }
 
  est.shd = lapply(1:(2*trials), est.shd, est.am.list, adjm.list)
  est.eqclass.shd =  lapply(1:(2*trials), eqclass.shd, est.am.list, adjm.list)
#   print("distance to ground truth estimated")
  
  
 
  stars.subs = lapply(subs.res, function(x) x[[1]])
  stars.resamp = lapply(resamp.res, function(x) x[[1]])
  
  subs.density = sapply(subs.res, function(x) x[[2]])
  resamp.density = sapply(resamp.res, function(x) x[[2]])
  
  true.nedges = sapply(adjm.list, function(x) sum(x))
  est.nedges = sapply(est.am.list, function(x) sum(x))
  
  export.data = data.frame(est.shd = as.numeric(est.shd),
                           est.eqclass.shd = as.numeric(est.eqclass.shd),
                           stars.subs = as.numeric(stars.subs),
                           stars.resamp = as.numeric(stars.resamp),
                           true.nedges= true.nedges,
                           est.nedges = est.nedges,
                           subs.nedges = choose(nnodes, 2)*subs.density,
                           resamp.density = choose(nnodes,2)*resamp.density,
                           data.size = effective.size.seq)
  
  write.table(export.data, 
              file = paste("est_data_nnodes_", nnodes,  "_",  method, "_nintervent_", nintervent, ".csv", sep = ""),
              row.names = F)
  
  cols = c(rep("green4", trials), rep("red", trials))
  #   plot(est.eqclass.dist, stars.subs,
  #        pch = 20, cex = 0.7, col = cols,
  #        ylab = "StARS",
  #        xlab = "sid eqclass distance",
  #        main = paste("stars subsample vs sid distance", method))
  #   plot(est.eqclass.diff, stars.subs, 
  #        pch = 20, cex = 0.7, col = cols,
  #        ylab = "StARS",
  #        xlab = "sid eqclass difference",
  #        main = paste("stars subsample vs sid difference", method))
  # 
  #   plot(est.eqclass.dist, stars.resamp,
  #        pch = 20, cex = 0.7, col = cols,
  #        ylab = "StARS",
  #        xlab = "sid eqclass distance",
  #        main = paste("stars resample vs sid distance", method))
  #   plot(est.eqclass.diff, stars.resamp,
  #        pch = 20, cex = 0.7,  col = cols,
  #        ylab = "StARS",
  #        xlab = "sid eqclass difference",
  #        main = paste("stars resample vs sid difference", method))
  plot(est.shd, stars.subs,
       pch = 20, cex = 0.7,  col = cols,
       ylab = "StARS",
       xlab = "shd to true graph",
       main = paste("stars subsample vs shd to true graph", method))
  plot(est.eqclass.shd, stars.subs,
       pch = 20, cex = 0.7,  col = cols,
       ylab = "StARS",
       xlab = "shd to true CPDAG",
       main = paste("stars subsample vs shd to true CPDAG", method))
  
  plot(est.shd, stars.resamp,
       pch = 20, cex = 0.7,  col = cols,
       ylab = "StARS",
       xlab = "shd to true graph",
       main = paste("stars resample vs shd to true graph", method))
  plot(est.eqclass.shd, stars.resamp,
       pch = 20, cex = 0.7,  col = cols,
       ylab = "StARS",
       xlab = "shd to true CPDAG",
       main = paste("stars resample vs shd to true CPDAG", method))
}


all.seq.sem = function(nnodes, densities, size.seq, effect.range, noise.range){
  
  sem.list = list()
  for(d in densities){
    sem.list = append(sem.list, replicate(length(size.seq), get.random.sem2(nnodes,d,effect.range, noise.range), simplify = F))
  }
  return(sem.list)
  
}

#all.seq.sem(densities, nnodes, size.seq, c(-3,3), c(0.1,2))

densities = seq(0.2,0.7,0.1)
nnodes = 11
nintervent = 7
size.seq = seq(30, 600, 15)
method = "gies"
dil.fact = 0.5

set.seed(6)
pdf(paste("stars_interventions_nnodes_", nnodes, "_nintervent_", nintervent, method, ".pdf", sep = ""))
score.intervention(nnodes, nintervent, 0, densities, size.seq, c(0.2,3), c(0.05,1.2), method, 200, dil.fact)
dev.off()

# 
# size.seq = c(50, 100, 150)
# nnodes = 8
# densities = c(0.3,0.4)
# 
# score.intervention(nnodes, 4, 0, densities, size.seq, c(-3,3), c(0.05, 1.2), "gies", 10, 0.5)
