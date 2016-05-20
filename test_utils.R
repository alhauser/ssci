################################################
# This scribt contains various functions used in
# several parts of my work
##############################################



fit.pcalg = function(data, method, method.name, indepTest = gaussCItest, alpha = 0.1, targets = NULL, target.index = NULL){
  
  # wraper for the estimations with the pcalg package
  # data: data frame
  # method: variable of estimation method in pcalg package
  # method.name: str of method name
  # return: the fitted object ( list of essgraph and repr)
  
  
  if(method.name %in% c("pc", "fci", "rfci")){
    suffStat = list(C = cor(data), n = nrow(data))
    fit = method(suffStat, indepTest, p = ncol(data), alpha)
  }else if(method.name == "ges"){
    score = new("GaussL0penObsScore", as.matrix(data))
    fit = method(ncol(data), score)
  } else if(method.name %in% c("gies", "simy")){
    int.score = new("GaussL0penIntScore", as.matrix(int.data), targets = targets, target.index = target.index)
    fit = method(ncol(data), targets, score)
  }
  return(fit)  
}
test.bn = function(data){
  #uses all the functions in bnlearn and prints the resulting graphs
  # data: data.frame
  
  bn.gs = gs(data)
  bn.iamb = iamb(data)
  bn.fast.iamb = fast.iamb(data)
  bn.inter.iamb = inter.iamb(data)
  bn.hc = hc(data)
  
  if(all.equal(bn.gs, bn.iamb) == T){
    print("gs equal iamb: True")
  }
  if(all.equal(bn.gs, bn.fast.iamb) == T){
    print("gs equal fast.iamb: True")
  }
  if(all.equal(bn.gs, bn.inter.iamb)== T){
    print("gs equal inter.iamb: True")
  }
  if(all.equal(bn.gs, bn.hc) != T){
    print("gs equal hc: False")
  }
  
  # printing graphs
  
  plot(bn.gs, main = "constrain based")
  plot(bn.hc, main = "hill climbing")
  
}
#-------------------------------------------------------------------------
#learning with pcalg and comparing all methods
#-------------------------------------------------------------------------

test.pcalg = function(data){
  #uses all the functions on observational data in pcalg and prints the resulting graphs
  # data: data.frame
  suffStat = list(C = cor(data), n = nrow(data))
  
  pc.fit = pc(suffStat, indepTest = gaussCItest, p = ncol(data), alpha = 0.01 )
  plot(pc.fit, main = "pc.fit")
  
  score = new("GaussL0penObsScore", as.matrix(data))
  
  ges.fit = ges(ncol(data), score)
  plot(ges.fit$essgraph, main = "ges.fit")
  
  fci.fit = fci(suffStat, indepTest = gaussCItest, p = ncol(data), alpha = 0.01 )
  plot(fci.fit, main = "fic.fit")
  
  rfci.fit = rfci(suffStat, indepTest = gaussCItest, p = ncol(data), alpha = 0.01 )
  plot(rfci.fit, main = "rfci.fit")
  
}

test.pcalg.int = function(int.data, int.targets, int.target.index){
  #uses all the functions on interventional data in pcalg and prints the resulting graphs
  # int.data : data frame
  # int.targets: list of intervention targets (number of variable in data frame
  #                                            0 for no intervetnion)
  # int.target.index: vector with index of intervention target in int.targets
  int.score = new("GaussL0penIntScore", as.matrix(int.data), targets = int.targets, target.index = int.target.index)
  
  gies.fit = gies(ncol(int.data), int.targets, int.score)
  plot(gies.fit$essgraph, main = "gies.fit")
  
  simy.fit = simy(ncol(int.data), int.targets, int.score)
  plot(simy.fit$essgraph, main = "simy.fit")
  
}

get.lin.gaus = function(sem = NULL, mmat = sem[[1]], sdev = sem[[2]], size){
  
  # mmat matrix with the coeficients bx_y
  # sdev sd of noise term
  # return : data.frame #size samples of the sem specified with mmat and sdev

  
  # initiate data frame with 0 for every variable in size rows
  data = data.frame(X1 = numeric(size))
  for(x in seq(ncol(mmat))){
    v = paste("X", x, sep = "")
    data[v] = rep(0, size)
  }
  
  #go through system by variables in causal order
  for(v in seq(nrow(mmat))){
    
    # calculate value of variable from previous ones, the n add
    # variables lower in causal order are still 0
    data[v] = as.matrix(data) %*% mmat[v,] + rnorm(size, 0, sdev[v]) 
  }
  return(data)
}

#----------------------------------------------------------------------------
# simulating interventions
#----------------------------------------------------------------------------

get.int.lin.gaus = function(mmat, target, tar.value, sdev, size){
  
  #mmat matrix with the coeficients bx_y
  #sdev sd of noise term
  # return : #size samples of the sem specified with mmat and sdev
  
  #initiate data frame with 0 for every variable in size rows
  data = data.frame(X1 = numeric(size))
  for(x in seq(ncol(mmat))){
    v = paste("X", x, sep = "")
    data[v] = rep(0, size)
  }
  
  # go through system by variables in causal order
  intervention = 1
  for(v in seq(ncol(data))){
    
    if(v %in% target){
      data[v] = rep(tar.value[intervention], size)
      intervention = intervention + 1
    }
    else{
      data[v] = as.matrix(data) %*% mmat[v,] + rnorm(size, 0, sdev[v]) 
    }  
  }
  return(data)
}

#-----------------------------------------------------------------------
# exponential noise
#-----------------------------------------------------------------------

get.lin.exp = function(mmat, lambda, size){
  
  #mmat matrix with the coeficients bx_y
  #lambda : data frame with upper and lower limit of noise for every node 
  # return : #size samples of the sem specified with mmat and sdev
  # set.seed(seed)
  
  #initiate data frame with 0 for every variable in size rows
  data = data.frame(X1 = numeric(size))
  for(x in seq(ncol(mmat))){
    v = paste("X", x, sep = "")
    data[v] = rep(0, size)
  }
  # go through system by variables in causal order
  for(v in seq(nrow(mmat))){
    data[v] = as.matrix(data) %*% mmat[v,] + rexp(size, lambda[v]) 
  }
  return(data)
}

get.int.exp = function(mmat, target, tar.value, lambda, size){
  
  #mmat matrix with the coeficients bx_y
  #sdev sd of noise term
  # return : #size samples of the sem specified with mmat and sdev
  
  #initiate data frame with 0 for every variable in size rows
  data = data.frame(X1 = numeric(size))
  for(x in seq(ncol(mmat))){
    v = paste("X", x, sep = "")
    data[v] = rep(0, size)
  }
  
  # go through system by variables in causal order
  intervention = 1
  for(v in seq(ncol(data))){
    
    if(v %in% target){
      data[v] = rep(tar.value[intervention], size)
      intervention = intervention + 1
    }
    else{
      data[v] = as.matrix(data) %*% mmat[v,] + rexp(size, lambda) 
    }  
  }
  return(data)
}

#-----------------------------------------------------------------------
# uniform noise
#-----------------------------------------------------------------------

get.lin.unif = function(mmat, rnge, size){
  
  #mmat matrix with the coeficients bx_y
  #rnge : data frame with upper and lower limit of noise for every node 
  # return : #size samples of the sem specified with mmat and sdev
  #set.seed(seed)
  

  #initiate data frame with 0 for every variable in size rows
  data = data.frame(X1 = numeric(size))
  for(x in seq(ncol(mmat))){
    v = paste("X", x, sep = "")
    data[v] = rep(0, size)
  }
  
  # go through system by variables in causal order
  for(v in seq(nrow(mmat))){
    data[v] = as.matrix(data) %*% mmat[v,] + runif(size, rnge[v,1], rnge[v,2]) 
  }
  return(data)
}

get.int.unif = function(mmat, target, tar.value, rnge, size){
  
  #mmat matrix with the coeficients bx_y
  #sdev sd of noise term
  # return : #size samples of the sem specified with mmat and sdev
  
  #initiate data frame with 0 for every variable in size rows
  data = data.frame(X1 = numeric(size))
  for(x in seq(ncol(mmat))){
    v = paste("X", x, sep = "")
    data[v] = rep(0, size)
    
  }
  
  intervention = 1
  for(v in seq(nrow(mmat))){
    
    if(v == target){
      data[v] = rep(tar.value[intervention], size)
      intervention = intervention + 1  
    }
    else{
      data[v] = as.matrix(data) %*% mmat[v,] + runif(size, rnge[v,1], rnge[v,2]) 
    }  
  }
  return(data)
}


#-----------------------------------------------------------------------
# non linear terms
#-----------------------------------------------------------------------

get.pot.gaus = function(mmat, sdev, pot, size){
  
  #mmat matrix with the coeficients bx_y
  #sdev sd of noise term
  # pot: matrix of exponents for the interactions 
  # return : #size samples of the sem specified with mmat and sdev
  
  #initiate data frame with 0 for every variable in size rows
  data = data.frame(X1 = numeric(size))
  for(x in seq(ncol(mmat))){
    v = paste("X", x, sep = "")
    data[v] = rep(0, size)
  }
  
  # calculate values
  for(i in seq(nrow(mmat))){
    #add effect defined by mmat and pot of every variable
    for(j in seq(ncol(mmat))){
      data[i] = data[i] + mmat[i,j] * data[j]^(pot[i,j])
    }
    # add noise
    data[i] = data[i] + rnorm(size, 0, sdev[i]) 
  }
  return(data)
  
}

get.int.pot = function(mmat, pot, target, tar.value, sdev, size){
  
  #mmat matrix with the coeficients bx_y
  #sdev sd of noise term
  # return : #size samples of the sem specified with mmat and sdev
  #set.seed(seed)
  
  #initiate data frame with 0 for every variable in size rows
  data = data.frame(X1 = numeric(size))
  for(x in seq(ncol(mmat))){
    v = paste("X", x, sep = "")
    data[v] = rep(0, size)
  }
  
  # calculate values
  intervention = 1
  for(i in seq(nrow(mmat))){
    
    if(i == target){
      data[i] = rep(tar.value[intervention], size)
      intervention = intervention + 1
    }
    else{
      #add effect defined by mmat and pot of every variable
      for(j in seq(ncol(mmat))){
        data[i] = data[i] + mmat[i,j] * data[j]^(pot[i,j])
      }
      # add noise
      data[i] = data[i] + rnorm(size, 0, sdev[i])   
    }
    
  }
  return(data)
  
}


get.lin.logis.gaus = function(mmat, logistic, l, k, sdev, size){
  #mmat matrix with the coeficients bx_y
  #sdev sd of noise term
  # return : #size samples of the sem specified with mmat and sdev
  # locistic: matrix , 1 if interaction is logistig, 0 otherwisw
  # l = vector: max of corresponding logistic interaction in topolocical order
  # k = vector: steepness of corresponding logistic interaction in topolocical order
  # sdev: vector standard variation of nois of each node
  # size: int number of values to simulate, sample size
  
  #initiate data frame with 0 for every variable in size rows
  data = data.frame(X1 = numeric(size))
  for(x in seq(ncol(mmat))){
    v = paste("X", x, sep = "")
    data[v] = rep(0, size)
  }
  
  # calculate values
  logistic.interaction = 1
  
  for(i in seq(nrow(mmat))){
    for(j in seq(ncol(mmat))){
      # calculate effect of j to i
      if(logistic [i,j] == 1){
        # logistic equation
        effect = l[logistic.interaction] / (1 + exp(-k[logistic.interaction] * data[j])) - 0.5 * l[logistic.interaction]
        logistic.interaction = logistic.interaction + 1
        #print(paste("logistic",i,j))
      }else{
        effect = mmat[i,j] * data[j]
        #print(paste("linear", i,j))
      }
      # add effect
      data[i] = data[i] + effect
    }
    #add noise
    data[i] = data[i] + rnorm(size, 0, sdev[i]) 
  }
  return(data)
  
}

get.int.logis = function(mmat, logistic, targets, target.values, l, k, sdev, size){
  #mmat matrix with the coeficients bx_y
  #sdev sd of noise term
  # return : #size samples of the sem specified with mmat and sdev
  # locistic: matrix , 1 if interaction is logistig, 0 otherwisw
  # l = vector: max of corresponding logistic interaction in topolocical order
  # k = vector: steepness of corresponding logistic interaction in topolocical order
  # sdev: vector standard variation of nois of each node
  # size: int number of values to simulate, sample size
  
  #initiate data frame with 0 for every variable in size rows
  data = data.frame(X1 = numeric(size))
  for(x in seq(ncol(mmat))){
    v = paste("X", x, sep = "")
    data[v] = rep(0, size)
  }
  
  # calculate values
  logistic.interaction = 1
  intervention = 1
  for(i in seq(nrow(mmat))){
    if(i %in% targets){
      data[i] = target.values[intervention]
      intervention = intervention +1
    }
    else{
      for(j in seq(ncol(mmat))){
        # calculate effect of j to i
        if(logistic [i,j] == 1){
          # logistic equation
          effect = l[logistic.interaction] / (1 + exp(-k[logistic.interaction] * data[j])) - 0.5 * l[logistic.interaction]
          logistic.interaction = logistic.interaction + 1
          #print(paste("logistic",i,j))
        }else{
          effect = mmat[i,j] * data[j]
          #print(paste("linear", i,j))
        }
        # add effect
        data[i] = data[i] + effect
      }
      #add noise
      data[i] = data[i] + rnorm(size, 0, sdev[i]) 
    }
  }
  return(data)
}
random.logistic.from.sem = function(sem, size, prob.logistic = 0.5, lim.range = c(6,14), kpos.range = c(0.5, 4)){
  #sem: list [[1]] model matrix of linear interactions [[2]] vector of standard variation of nois of every node
  #size: int number of values to simulate, sample size
  # prob.logistic: float (0,1) probability for a interaction to be logistig
  # lim.range: float (0,inf) range for limit parameter in logistic equation
  # kpos.range: float (0, inf) positive range for k (slope in logistic eq) values will be multiplied by -1 with prob 0.5
  
  mmat = sem[[1]]
  sdev = sem[[2]]
  
  # estimate which interactions should be 
  interactions = sum(mmat != 0)
  log.interactios = rbinom(interactions, 1, prob.logistic)
  log.mat = mmat
  log.mat[which(log.mat!= 0)] = log.interactios
  
  # get parameters for logistic interactions
  nlog = sum(log.interactios)
  lim = runif(nlog, lim.range[1], lim.range[2])
  k = runif(nlog, kpos.range[1], kpos.range[2])
  neg = rbinom(nlog, 1, 0.5)
  k[which(neg == 1)] = -k[which(neg == 1)]
  
  return(get.lin.logis.gaus(mmat, log.mat, lim, k, sdev, size))
}

random.logistic.list = function(index, sem.list, size.list, prob.logistic = 0.5, lim.range = c(6,14), kpos.range = c(0.5, 4)){
  # wrapper for applying random.logistic.from.sem on a list of sem and list of sizes
  # index: int in 1:length sem.list
  # sem.list: list of sem
  # size. list: list of corresponding sizes
  return(random.logistic.from.sem(sem.list[[index]], size.list[[index]], prob.logistic, lim.range, kpos.range))  
}

get.random.sem = function(n, density, effect.rnge, noise.rnge){
  # generates a random sem
  # n: int, number of nodes
  # density: probability for a possible edge to be non 0
  # effect.rnge: lower and upper limit for edge weight
  # noise.rnge: lower and upper limit for noise standard deviation
  # return : list of 2 mmat defining the interactions, sdev standart deviation of noises
  
  # estimate edge values
  poss.edges = (n-1)*n/2 # chose(n, 2) lower triangle of matrix
  realised.edges = rbinom(poss.edges, 1, density)# indicator for edges to set to zero 
  edges = runif(poss.edges, effect.rnge[1], effect.rnge[2] )
  edges = edges * realised.edges # set selcted edges to 0
  
  mmat = matrix(0,n,n)
  c = 1
  # write edges into mmat (out going [i,j] is edge form i to j)
  for(i in 1:n){
    for(j in 1:n){
      if(i > j){
        mmat[i,j] = edges[c]
        c = c + 1
      }
    }
  }
  sdev = runif(n, noise.rnge[1], noise.rnge[2])
  return(list(mmat, sdev))
}


get.random.sem2 = function(n, density, pos.effect.rnge, noise.rnge){
  # generates a random sem
  # n: int, number of nodes
  # density: probability for a possible edge to be non 0
  # pos.effect.rnge: lower and upper limit for edge weight
  #                  half of the edgeweights will be multipied by -1 
  #                  this is to avoid edgeweights close to 0
  # noise.rnge: lower and upper limit for noise standard deviation
  # return : list of 2 mmat defining the interactions, sdev standart deviation of noises
  
  # estimate edge vauels
  poss.edges = (n-1)*n/2
  edges = runif(poss.edges, pos.effect.rnge[1], pos.effect.rnge[2]) * rbinom(poss.edges, 1, density)
  # make edges negative with prob 0.5
  nedges = sum(edges > 0)
  neg.edges = rbinom(sum(edges > 0), 1, 0.5)
  neg.edges[which(neg.edges == 0)] = -1
  edges[edges != 0] = edges[edges != 0] * neg.edges
  
  mmat = matrix(0,n,n)
  c = 1
  #write edges into mmat
  for(i in 1:n){
    for(j in 1:n){
      if(i > j){
        mmat[i,j] = edges[c]
        c = c + 1
      }
    }
  }
  sdev = runif(n, noise.rnge[1], noise.rnge[2])
  return(list(mmat, sdev))
}
adjm.from.mmat = function(mmat){
  
  true.am = t(mmat)#mmat contains outgoing edges, incomming ones are needed used for adjiacency matrixes in R
  true.am[which(true.am != 0)] = 1
  
  # add nodenames (needed for get.adjm.from.bn)
  colnames(true.am) = rownames(true.am) = paste(rep("X", dim(true.am)[1]), seq(dim(true.am)[1]),sep = "")
  # named X1, X2, ..., Xn
  return(true.am)
}

graph.from.mmat = function(mmat){
  # mmat: model matrix
  # return: igraph object of model
  g = graph.adjacency(adjm.from.mmat(mmat))
}


get.adjm.from.essgraph = function(essgraph){
  
  # essgraph: essgraph objecto of pcalg fit
  # return: adjacency matrix of essgraph
  
  nnodes = length(essgraph$.in.edges)
  am = matrix(rep(0, nnodes^2), nrow = nnodes)
  c= 1
  for( i in essgraph$.in.edges){
    for(r in i){
      am[r,c] = 1
    }
    c = c + 1
  }
  return(am)
}

get.adjm.from.bn = function(fit, cpdag = F)
  # bn objects do not contain a graph object ! =((
  # returns a adjancy matrix of all variables including unconnected ones
{
  tmp.am = get.adjacency(graph.edgelist(as.matrix(fit$arcs), directed = T))
  # does not include nodes with no conections
  
  am = matrix(0, nnodes(fit), nnodes(fit))
  colnames(am) = rownames(am) = nodes(fit)
  
  # am with only 0, same row and columnames as true.am
  # write all values from tmp.am to corresponding entry in am (same row and collumnames)
  for(i in rownames(am)){
    for(j in colnames(am)){
      if(i %in% rownames(tmp.am) && j %in% colnames(tmp.am)){
        am[i,j] = tmp.am[i,j]
      }
    }
    
  }
  
  if(cpdag){
    g = graphAM(am, edgemode = "directed")
    g = as(g, "graphNEL")
    am = dag2cpdag(g)
    am = as(am , "matrix")
  }
  return(am)
}


scramble.dat.list = function(index, df.list, perm.list){
  # index: int in 1:length(df.list.index)
  # df.list: lits of data frames
  # perm.list: list of permutations of ncol(df.list)
  # return: columne wise permutated data frame of df.list[[index]]
  return(df.list[[index]][perm.list[[index]]])
} 

scramble.adjmat.list = function(index, perm.list, adjm.list){
  # index: int in 1:length(df.list.index)
  # adjm.list: lits of data adjacency matrixes
  # perm.list: list of permutations of ncol(adjm.list)
  # return: permutated adjacancy matrix of adjm.list[[index]]
  perm = perm.list[[index]]
  adjm = adjm.list[[index]]
  return(adjm[perm, perm])
}
scramble.data= function(df, perm = sample(1:ncol(df), replace = F)){
  return(df[perm])
}

true.eq.upper = function(true.am){
  # true.am: adjacency matrix of a true DAG
  # return: upper limit of sid of a graph in the true cpdag/eqclass to true graph
  if(length(colnames(true.am)) == 0){
    colnames(true.am) = 1:ncol(true.am)
  }
  # get cpdag
  g = graphAM(true.am, edgemode = "directed")
  g = as(g, "graphNEL")
  am = dag2cpdag(g)
  
  return(structIntervDist(true.am, am)$sidUpperBound)
}

cpdag.from.adjm = function(adjm, matrix = T){
  # adjm: adjacency matrix
  # matrix: boolean if TRUE only adjacency matrix is returned
  # return: matrix or graphAM object of corresponding cpdag
  if(length(colnames(adjm)) == 0){
    colnames(adjm) = 1:ncol(adjm)
  }
  g = graphAM(adjm, edgemode = "directed")
  g = as(g, "graphNEL")
  cpdag = as(dag2cpdag(g), "graphAM")
  if(matrix){
    cpdag = cpdag@adjMat
  }
  return(cpdag)
}

exp.from.list =function(sem, size){ 
  # sem: list of two, matrix of modle, vector of standard deviations
  # size: sample size 
  # return: data.frame with size samples from sem
  mmat = sem[[1]]
  lambda = 1/sem[[2]]
  d = get.lin.exp(mmat, lambda, size)
  return(d)
}

exp.random.size =function(sem, upper, lower){
  # sem: list of two, matrix of modle, vector of standard deviations
  # size: upper: int upper limit for sample size
  # size: lower: int lower limit for sample size
  # return: data frame 
  
  mmat = sem[[1]]
  lambda = 1/sem[[2]]
  size = runif(1, lower, upper)
  d = get.lin.exp(mmat, lambda, size)
  return(d)
}

exp.size.seq = function(index, sem.list, size.seq){
  # index: int in 1:length(sem.list)
  # sem.list: list of sem
  # size.seq : list of int sample size
  # return: data frame with size.seq[[index]] samples of sem[[index]]
  mmat = sem.list[[index]][[1]]
  lambda = 1/sem.list[[index]][[2]]
  size = size.seq[[index]]
  return(get.lin.exp(mmat, lambda, size))
}

lin.gaus.random.size = function(sem, upper, lower){
  
  # sem: list of two, matrix of modle, vector of standard deviations
  # size: upper: int upper limit for sample size
  # size: lower: int lower limit for sample size
  # return: data frame 
  size = runif(1, lower, upper)
  return(get.lin.gaus(sem, size = size))  
}

lin.gaus.size.seq = function(index, sem.list, size.seq ){
  # index: int in 1:length(sem.list)
  # sem.list: list of sem
  # size.seq : list of int sample size
  # return: data frame with size.seq[[index]] samples of sem[[index]]
  return(get.lin.gaus(sem = sem.list[[index]], size = size.seq[index]))
}

est.adjm = function(sample, method, environement = NULL, interventions = NULL, cpdag = NULL){
  #sample: data frame of data to be fitted
  #method: str method name of a algorithm in bnlearn of CompareCausalNetworks
  # environement: vector of nrow(sample), see getParents 
  # interventions : list of nrow(sampel), see getParents
  #return: estimated adjacency matrix 
  
  # methods in CompareCausalNetworks
  if(method %in% c("ICP", "hiddenICP", "backShift", "pc","LINGAM", "ges", "gies",
                   "CAM", "rfci", "regression", "bivariateANM","bivariateCAM")){
    am = as.matrix(getParents(sample, method = method, directed = F, environement, interventions))  
  }
  # methods in bnlearn
  else{
    fit = switch(method,
                 hc = {hc(sample)},
                 gs = {gs(sample)},
                 iamb = {iamb(sample)},
                 fast.iamb = {fast.iamb(sample)},
                 intrer.iamb = {inter.iamb(sample)},
                 hc = {hc(sample)},
                 tabu = {tabu(sample)},
                 mmhc = {mmhc(sample)},
                 rsmax2 = {rsmax2(sample)},
                 stop("No valid estimation method given. Method must be implemented in the package bnlearn or CompareCausalNetworks!")
    )
    am = get.adjm.from.bn(fit, cpdag = (method == "hc"))
  }
  return(am)
}