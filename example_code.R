
source("test_utils.R")
source("scoring.R")
set.seed(6)
#generate a random sem for simulations with get.random.sem2 @ test_utils.R
sem = get.random.sem2(n = 7, density = 0.3, pos.effect.rnge = c(0.2,3), noise.rnge = c(0.1, 1.2))

#get graph of generated sem
true.am = adjm.from.mmat(sem[[1]]) # @ test_utils.R
true.graph = graphAM(true.am, edgemode = "directed")

#create a plot with sid vs sample size with linear gaussian data and ges algorithm
size = seq(50, 2000, 20)
dat.list = lapply(size, FUN = get.lin.gaus, mmat = sem[[1]], sdev = sem[[2]], sem = NULL)
source("sid/sid_samplesize.R")
plot.sid(dat.list, size, ges, "ges", true.am )

# compare estimates on linear gaussian data to estimates on linear exponential data
par(mfrow = c(1,2))
plot(true.graph)
gauss.data = get.lin.gaus(sem = sem, size = 1000)
est.am = getParents(gauss.data, method = "ges", directed = F)
plot(graphAM(as.matrix(est.am), edgemode = "directed"))
est.stars.score(gauss.data, estimation.method = "ges") # @ scoring.R

exp.data = get.lin.exp(mmat = sem[[1]], lambda = 1/sem[[2]], size = 1000)
est.am = getParents(exp.data, method = "ges", directed = F)
plot(true.graph)
plot(graphAM(as.matrix(est.am), edgemode = "directed"))
est.stars.score(dat= exp.data, estimation.method = "ges")

experiment.list = list()
experiment.list[[1]] = list(get.lin.gaus(sem = sem, size= 300), integer(0))
experiment.list[[2]] = list(get.int.lin.gaus(mmat= sem[[1]], target = 2, tar.value = 0, sdev = sem[[2]], size = 300),2) 
ful.dat = concatenate.experiments(experiment.list)

getParents(as.data.frame(ful.dat[[1]]), method = "gies", interventions = ful.dat[[2]], directed = F)
