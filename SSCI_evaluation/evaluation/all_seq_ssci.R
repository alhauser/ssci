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
  shd = sum(est.adjm.list[[index]] != cpdag)
  return(shd)
}

est.sid.measures = function(index, est.adjm.list, true.adjm.list){
  sid = structIntervDist(true.adjm.list[[index]], est.adjm.list[[index]])
  true.upper = true.eq.upper(true.adjm.list[[index]])
  dist = sid$sidLowerBound + abs(sid$sidUpperBound - true.upper)
  diff = sid$sidLowerBound + sid$sidUpperBound - true.upper
  return(list(dist, diff))
}


score.all.seq =  function(nnodes, densities, size.seq, effect.rnge, noise.rnge, method, resamps, dil.fact){
  #n: int number of variables in sem, nodes in graph
  #density: R[0,1] density of graph, fraction of edges to number of possible edges
  #effect.rnge: vector of upper and lower limit to sample magnitude of causal effects from
  #noise.rnge: vector of upper and lower limit oto sample the standardvariation of the noise terms
  #trails: int numer of different sem to test
  #method: str name of algo to infere model with ( CompareCausalNetworks or bnlearn)
  #resamp.method: str "subsample" or "resample"
  #resamps: int number of resamples used for score estimation
  #size: int number of samples drawn from every sem
  

  sem = all.seq.sem(nnodes, densities, size.seq, effect.rnge, noise.rnge)
  sem = append(sem, all.seq.sem(nnodes, densities, size.seq, effect.rnge, noise.rnge))
  mmat.list = lapply(sem, function(x) x[[1]])
  #sdev = lapply(sem, function(x) x[[2]])
  adjm.list = lapply(mmat.list, adjm.from.mmat)
  
#   perm.list = replicate(trials, sample(1:n, replace = F), simplify = F)
#   adjm.list = lapply(1:trials, scramble.adjmat.list, perm.list, adjm.list)
#   
  
  effective.size.seq = rep(rep(size.seq, length(nnodes)*length(densities)),2)
  gaus.index = 1:(length(nnodes)*length(densities)*length(size.seq))
  exp.index = (length(gaus.index)+1):(length(effective.size.seq))
  trials = length(effective.size.seq)

  lin.gaus.data = lapply(gaus.index, lin.gaus.size.seq,  sem, effective.size.seq)
  #lin.gaus.data = lapply(gaus.index, scramble.dat.list, lin.gaus.data, perm.list)
  
  lin.exp.data = lapply(exp.index, exp.size.seq, sem, effective.size.seq)
  #exp.perm.list = perm.list[exp.index]
  #lin.exp.data = lapply(1:length(lin.exp.data), scramble.dat.list, lin.exp.data, exp.perm.list)
  dat = append(lin.gaus.data, lin.exp.data)
  #dat = lapply(1:trials, scramble.dat.list, dat, perm.list)
  
#   random.dags = replicate(100, randomDAG(n, density))
#   score.random = stars.score(random.dags)
#   
  est.am.list = lapply(dat, est.adjm, method = method)
  
#   sid.measures = lapply(1:trials, est.sid.measures, est.am.list, adjm.list)
#   
#   est.eqclass.dist = lapply(sid.measures, function(x) x[[1]])
#   est.eqclass.diff = lapply(sid.measures, function(x) x[[2]])
# 
#   est.eqclass.dist = lapply(1:trials, est.eqclass.dist, dat, adjm.list, method)
#   est.eqclass.diff = lapply(1:trials, est.eqclass.diff, dat, adjm.list, method)
#   

  est.shd = lapply(1:trials, est.shd, est.am.list, adjm.list)
  est.eqclass.shd =  lapply(1:trials, eqclass.shd, est.am.list, adjm.list)
  print("distance to ground truth estimated")
  
  
  subs.res = lapply(dat, score.edgedist, "subsample", method, resamps, dil.fact = dil.fact)
  print("stars subsampling estimated")
  #resamp.res = lapply(dat, score.edgedist, "resample", method, resamps)
  print("stars resampling estimated")
  stars.subs = lapply(subs.res, function(x) x[[1]])
  #stars.resamp = lapply(resamp.res, function(x) x[[1]])
  
  subs.density = sapply(subs.res, function(x) x[[2]])
  #resamp.density = sapply(resamp.res, function(x) x[[2]])
  
  true.nedges = sapply(adjm.list, function(x) sum(x))
  est.nedges = sapply(est.am.list, function(x) sum(x))
  
  export.data = data.frame(est.shd = as.numeric(est.shd),
                           est.eqclass.shd = as.numeric(est.eqclass.shd),
                           stars.subs = as.numeric(stars.subs),
                          # stars.resamp = as.numeric(stars.resamp),
                           true.nedges= true.nedges,
                           est.nedges = est.nedges,
                           subs.nedges = choose(nnodes, 2)*subs.density,
                           #resamp.density = choose(nnodes,2)*resamp.density,
                           data.size = effective.size.seq)
  
  write.table(export.data, file = paste("est_data_nnodes_", method, "_", nnodes, ".txt", sep = ""), row.names = F)
  
  cols = c(rep("green4", length(lin.gaus.data)), rep("red", length(lin.exp.data)))
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

#   plot(est.shd, stars.resamp,
#        pch = 20, cex = 0.7,  col = cols,
#        ylab = "StARS",
#        xlab = "shd to true graph",
#        main = paste("stars resample vs shd to true graph", method))
#   plot(est.eqclass.shd, stars.resamp,
#        pch = 20, cex = 0.7,  col = cols,
#        ylab = "StARS",
#        xlab = "shd to true CPDAG",
#        main = paste("stars resample vs shd to true CPDAG", method))
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
nnodes = 15 
size.seq = seq(100, 1500, 50)
method = "ges"
dil.fact = 0.5

set.seed(6)
pdf(paste("stars_all_seq_nnodes_", nnodes, "_", method, ".pdf", sep = ""))
score.all.seq(nnodes, densities, size.seq, c(0.2,3), c(0.05,1.2), method, 200, dil.fact)
dev.off()


# size.seq = c(50, 100, 150)
# nnodes = 4
# densities = c(0.3,0.4)
# 
# 
# score.all.seq(nnodes,densities, size.seq, c(0.2,3), c(0.05,1.2), "ges", 10, 0.5)
