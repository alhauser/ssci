source("../../test_utils.R")
source("../sid_samplesize.R")
set.seed(42)#setting the same seed will result in very simmilar structures in the SEM

sem = get.random.sem2(12, 0.3, c(0.3,3), c(0.1,1.2)) # small but dense

mmat = sem[[1]]
sdev = sem[[2]]

write.table(mmat, "mmat.txt",row.names = F, col.names = F)
write.table(sdev, "sdev.txt", row.names = F, col.names = F)
# 
# pdf("distribution.pdf")
# pairs(get.lin.gaus(mmat, sdev, 2000))
# dev.off()

#---------------------------------------------------------------
# create true adjacency matrix from mmat
#---------------------------------------------------------------
true.am= t(sem[[1]])#mmat is contains outgoing edges, incomming ones are needed for sid
true.am[which(true.am != 0)] = 1

# add nodenames (needed for get.adjm.from.bn)
colnames(true.am) = rownames(true.am) = paste(rep("X", dim(true.am)[1]), seq(dim(true.am)[1]), sep = "")
# named X1, X2, ..., Xn

g = graph.adjacency(true.am)
pdf("graph.pdf")
plot(g)
dev.off()

get.int.index = function(size, targets){
  
  obs = length(targets) + 1
  int.target.index = rep(1:obs, rep(size, obs))
  
}

size = seq(50, 4000, 10)

targets = c(3,5,6,10)
tar.values = c(1,1,1,1)
int.targets = list(integer(0), 3, 5, 6, 10)
int.size = round(size/(length(targets) +1))
int.index = lapply(int.size, get.int.index, targets)

source("../sid_lin_gauss.R")
source("../sid_unif_noise.R")
source("../sid_exp_noise.R")

quad.density = 0.3
quad = t(true.am)
quad[which(quad == 1)] = 1 + rbinom(sum(quad), 1, quad.density)

qubic = quad
qubic[which(qubic == 2)] = 3

source("../sid_non_lin.R")


write.table(quad, "quad.txt", col.names = F, row.names = F)
write.table(qubic, "cubic.txt", col.names = F, row.names = F)

logistic = quad
logistic[which(logistic != 0)] = logistic[which(logistic != 0)]-1

set.seed(6)
l = runif(sum(logistic), 3,8)
k = runif(sum(logistic), 0.1, 1.7)
neg = rbinom(sum(logistic), 1, 0.5)
k[which(neg == 1)] = -k[which(neg == 1)]

write.table(logistic, "logistic.txt", col.names = F, row.names = F)
write.table(l, "logis-l.txt", col.names = F, row.names = F)
write.table(k, "logis-k.txt",col.names = F, row.names = F)

source("../sid_logistic.R")
source("interactions.R")

