
source("../sid_samplesize.R")


get.interaction.data = function(size, mmat, sdev){
  
  data = data.frame(X1 = numeric(size))
  for(x in seq(ncol(mmat))){
    v = paste("X", x, sep = "")
    data[v] = rep(0, size)
  }
  
  data$X1 = rnorm(size, 0, sdev[1])
  data$X2 = rnorm(size, 0, sdev[2])
  data$X3 = mmat[3,2] * data$X2 + rnorm(size, 0, sdev[3])
  data$X4 = mmat[4,1] * data$X1 * data$X3 + rnorm(size, 0, sdev[4]) #interaction
  data$X5 = mmat[5,3] * data$X3 * data$X4 + rnorm(size, 0, sdev[5]) #interaction
  data$X6 = mmat[6,4] * data$X4 + mmat[6,5] * data$X5 + rnorm(size, 0, sdev[6])
  data$X7 = mmat[7,1] * data$X1 * data$X4 + mmat[7,6] * data$X6 + rnorm(size, 0, sdev[7]) # interaction
  
  #3 interaction terms
  
  return(data)
  
}
dat.list = lapply(size, FUN = get.interaction.data, mmat = mmat, sdev = sdev)

get.int4 = function(size, mmat, sdev){
  
  data = data.frame(X1 = numeric(size))
  for(x in seq(ncol(mmat))){
    v = paste("X", x, sep = "")
    data[v] = rep(0, size)
  }
  
  data$X1 = rnorm(size, 0, sdev[1])
  data$X2 = rnorm(size, 0, sdev[2])
  data$X3 = mmat[3,2] * data$X2 + rnorm(size, 0, sdev[3])
  data$X4 = 0
  data$X5 = mmat[5,3] * data$X3 * data$X4 + rnorm(size, 0, sdev[5])
  data$X6 = mmat[6,4] * data$X4 + mmat[6,5] * data$X5 + rnorm(size, 0, sdev[6])
  data$X7 = mmat[7,1] * data$X1 * data$X4 + mmat[7,6] * data$X6 + rnorm(size, 0, sdev[7])
  
  #3 interaction terms
  
  return(data)
}


get.int3= function(size, mmat, sdev){
  
  data = data.frame(X1 = numeric(size))
  for(x in seq(ncol(mmat))){
    v = paste("X", x, sep = "")
    data[v] = rep(0, size)
  }
  
  data$X1 = rnorm(size, 0, sdev[1])
  data$X2 = rnorm(size, 0, sdev[2])
  data$X3 = 0
  data$X4 = mmat[4,1] * data$X1 * data$X3 + rnorm(size, 0, sdev[4])
  data$X5 = mmat[5,3] * data$X3 * data$X4 + rnorm(size, 0, sdev[5])
  data$X6 = mmat[6,4] * data$X4 + mmat[6,5] * data$X5 + rnorm(size, 0, sdev[6])
  data$X7 = mmat[7,1] * data$X1 * data$X4 + mmat[7,6] * data$X6 + rnorm(size, 0, sdev[7])
  
  #3 interaction terms
  
  return(data)
}
get.int.data = function(mmat, sdev, size){
  
  data = get.interaction.data(size, mmat, sdev)
  
  data = rbind(data, 
               get.int3(size, mmat, sdev),
               get.int4(size, mmat, sdev)) 
  
  return(data)
}

int.data = lapply(int.size, get.int.data, mmat = mmat, sdev = sdev)


pdf("sid_interactions.pdf")
#plot.sid(dat.list, size, hc, "hc", mmat)
#plot.sid(dat.list, size, gs, "gs", mmat)   
print("start pc")
plot.sid(dat.list, size, pc, "pc", true.am)
#plot.sid(dat.list, size, fci, "fci", mmat)
print("start ges")
plot.sid(dat.list, size, ges, "ges", true.am)
#plot.sid(dat.list, size, iamb, "iamb", mmat)
print("start gies")
plot.int.sid(int.data, int.index, int.targets, method = gies, method.name = "gies", true.am, int.size)

dev.off()

