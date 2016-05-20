source("../../scoring.R")
source("../../test_utils.R")
library(ggplot2)
set.seed(6)

random.dags.score = function(n, density, reps = 1000){
  # n: int number of nodes in DAG
  # density: float [0,1], density of DAG
  # reps : int, number of replicats
  # return: estimated SSCI of random Dags
  random.dags = replicate(reps, randomDAG(n, density))
  score = stars.score(random.dags)
  return(score)
} 

random.cpdag.score = function(n, density, reps){
  # n: int number of nodes in DAG
  # density: float [0,1], density of DAG
  # reps : int, number of replicats
  # return: estimated SSCI of random CPDAGs
  random.dags = replicate(reps, as.matrix(randomDAG(n, density)), simplify = F)
  random.cpdags = lapply(random.dags, cpdag.from.adjm)
  
  score = stars.score(random.cpdags)
  return(score)
  
}


density.plot = function(density = seq(0.1, 1, by = 0.06), nnodes, ndags, replicates){
# creates plots with SSCI of random DAGs as function of the density of the DAGs   
  
  df = data.frame(density = density)
  
  for(i in 1:length(density)){
    scores = replicate(replicates, random.dags.score( n = nnodes, density = density[i], reps = ndags))
    df$mean[i] = mean(scores)
    df$sdev[i] = sd(scores)
  }
  
  ggplot(df, aes(x = density, y = mean)) +ylab("SSCI") + 
    geom_errorbar(aes(ymin =mean - sdev, ymax = mean + sdev), width = 0.03) +
    geom_line() + geom_point()+
    ggtitle(paste("DAG nnodes = ", nnodes, ", ", "ndags = ", ndags, ", ", "replicates = ", replicates, sep = ""))
} 

cpdag.density.plot = function(density = seq(0.1, 1, by = 0.06), nnodes, ndags, replicates){
  # creates plots with SSCI of random DAGs as function of the density of the DAGs   
  
  df = data.frame(density = density)
  
  for(i in 1:length(density)){
    scores = replicate(replicates, random.cpdag.score( n = nnodes, density = density[i], reps = ndags))
    df$mean[i] = mean(scores)
    df$sdev[i] = sd(scores)
  }
  
  ggplot(df, aes(x = density, y = mean)) +ylab("SSCI") + 
    geom_errorbar(aes(ymin =mean - sdev, ymax = mean + sdev), width = 0.03) +
    geom_line() + geom_point()+
    ggtitle(paste("CPDAG nnodes = ", nnodes, ", ", "ncpdags = ", ndags, ", ", "replicates = ", replicates, sep = ""))
} 

pdf("density_plots.pdf")
#density.plot(nnodes = 5, ndags = 100, replicates = 5)
#density.plot(nnodes = 10, ndags = 100, replicates = 5)
density.plot(nnodes = 15, ndags = 200, replicates = 10)
cpdag.density.plot(nnodes = 15, ndags = 200, replicates = 10)
dev.off()


ndags.plot = function(ndags = seq(2, 400, by = 40), nnodes , density, replicates){
# creates plots with the SSCI of random dags as the function of the number of DAGs used for estimation
  
  df = data.frame(ndags = ndags)
  
  for(i in 1:length(ndags)){
    scores = replicate(replicates, random.dags.score( n = nnodes, density = density, reps = ndags[i]))
    df$mean[i] = mean(scores)
    df$sdev[i] = sd(scores)
  }
  
  ggplot(df, aes(x = ndags, y = mean)) +ylab("SSCI") + xlab("number of DAGs") +
    geom_errorbar(aes(ymin =mean - sdev, ymax = mean + sdev), width = 0.03) +
    geom_line() + geom_point() +
    ggtitle(paste("DAG nnodes = ", nnodes, ", ", "density = ", density, ", ", "replicates = ", replicates, sep = ""))
  
}


cpdag.ndags.plot = function(ncpdags = seq(2, 400, by = 40), nnodes , density, replicates){
  # creates plots with the SSCI of random cpdags as the function of the number of CPDAGs used for estimation
  
  df = data.frame(ncpdags = ncpdags)
  
  for(i in 1:length(ncpdags)){
    scores = replicate(replicates, random.cpdag.score(n = nnodes, density = density, reps = ncpdags[i]))
    df$mean[i] = mean(scores)
    df$sdev[i] = sd(scores)
  }
  
  ggplot(df, aes(x = ncpdags, y = mean)) +ylab("SSCI") + xlab("number of CPDAGs") +
    geom_errorbar(aes(ymin = mean - sdev, ymax = mean + sdev), width = 0.03) +
    geom_line() + geom_point() +
    ggtitle(paste("CPDAG nnodes = ", nnodes, ", ", "density = ", density, ", ", "replicates = ", replicates, sep = ""))
  
}

pdf("ncpdag25.pdf")
#ndags.plot(nnodes = 5, density = 0.3, replicates = 10)
#ndags.plot(nnodes = 10, density = 0.3, replicates = 10)
ndags.plot(nnodes = 7, density = 0.3, replicates = 10)
cpdag.ndags.plot(nnodes = 7, density = 0.3, replicates = 10)
#ndags.plot(nnodes = 25, density = 0.3, replicates = 10)
dev.off()

ndags.plot(ndags = seq(2, 500, by = 100), nnodes = 5, density = 0.1, replicates = 5 )
cpdag.ndags.plot(ncpdags = seq(2, 500, by = 100), nnodes = 5, density = 0.1, replicates = 5 )


nnodes.plot = function(ndags = 200, nnodes =  seq(3, 30, by = 3), density, replicates){
# creates plots with the SSCI as function of the number of nodes in the DAG
  
  df = data.frame(nnodes = nnodes)
  
  for(i in 1:length(nnodes)){
    scores = replicate(replicates, random.dags.score( n = nnodes[i], density = density, reps = ndags))
    df$mean[i] = mean(scores)
    df$sdev[i] = sd(scores)
  }
  
  ggplot(df, aes(x = nnodes, y = mean)) +ylab("SSCI") + xlab("number of nodes") +
    geom_errorbar(aes(ymin =mean - sdev, ymax = mean + sdev), width = 0.03) +
    geom_line() + geom_point() +
    ggtitle(paste("DAG ndags = " , ndags, ", density = ", density, ", replicates = ", replicates, sep = ""))
  
}

cpdag.nnodes.plot = function(ncpdags = 200, nnodes =  seq(3, 30, by = 3), density, replicates){
  # creates plots with the SSCI as function of the number of nodes in the CPDAG
  
  df = data.frame(nnodes = nnodes)
  
  for(i in 1:length(nnodes)){
    scores = replicate(replicates, random.cpdag.score( n = nnodes[i], density = density, reps = ncpdags))
    df$mean[i] = mean(scores)
    df$sdev[i] = sd(scores)
  }
  
  ggplot(df, aes(x = nnodes, y = mean)) +ylab("SSCI") + xlab("number of nodes") +
    geom_errorbar(aes(ymin =mean - sdev, ymax = mean + sdev), width = 0.03) +
    geom_line() + geom_point() +
    ggtitle(paste("CPDAG ncpdags = ", ncpdags, ", density = ", density, ", replicates = ", replicates, sep = ""))
  
}

nnodes.plot(ndags = 20, nnodes = c(5,8), density = 0.3, replicates = 5)
cpdag.nnodes.plot(ncpdags = 20, nnodes = c(5,8), density = 0.3, replicates = 5 )
pdf("nnodes_plots.pdf")
nnodes.plot(density = 0.3, replicates = 10)
#nnodes.plot(density = 0.2, replicates = 10)
cpdag.nnodes.plot(density = 0.3, replicates = 10)
#cpdag.nnodes.plot(density = 0.2, replicates = 10)
dev.off()
