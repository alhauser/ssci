
source("../sid_samplesize.R")

set.seed(6)
get.interaction.data = function(size, mmat, sdev){
  
  data = data.frame(X1 = numeric(size))
  for(x in seq(ncol(mmat))){
    v = paste("X", x, sep = "")
    data[v] = rep(0, size)
  }
  
  data$X1 = rnorm(size, 0, sdev[1])
  data$X2 = rnorm(size, 0, sdev[2])
  data$X3 = mmat[3,1] * data$X1 + rnorm(size, 0, sdev[3]) 
  data$X4 = rnorm(size, 0, sdev[4])
  data$X5 = mmat[5,4] * data$X4 + rnorm(size, 0, sdev[5])
  data$X6 = rnorm(size, 0, sdev[6]) 
  data$X7 = mmat[7,4] * data$X4 +  rnorm(size, 0, sdev[7])
  data$X8 = mmat[8,7] *  data$X7 + rnorm(size, 0, sdev[8]) 
  data$X9 = mmat[9,1] * data$X1 * data$X2 + mmat[9,5] *data$X5 +rnorm(size, 0, sdev[9])#interaction
  data$X10 = mmat[10,3]*data$X3 * data$X4 + mmat[10,5] * data$X5 * data$X6 + rnorm(size, 0, sdev[10])#interaction
  data$X11 = mmat[11,1] * data$X1 + mmat[11,2] * data$X2 * data$X9 + rnorm(size, 0, sdev[11]) # interaction
  data$X12= mmat[12,4] * data$X4 + mmat[12,6] * data$X6 * data$X7 + mmat[12,9] * data$X9 +
            mmat[12,10] * data$X10 + rnorm(size, 0, sdev[12]) # interaction
  
  #5 interaction terms
  
  return(data)
  
}

dat.list = lapply(size, FUN = get.interaction.data, mmat = mmat, sdev = sdev)


get.int3= function(size, mmat, sdev){
  
  
  data = data.frame(X1 = numeric(size))
  for(x in seq(ncol(mmat))){
    v = paste("X", x, sep = "")
    data[v] = rep(0, size)
  }  
  
  data$X1 = rnorm(size, 0, sdev[1])
  data$X2 = rnorm(size, 0, sdev[2])
  data$X3 = 1 
  data$X4 = rnorm(size, 0, sdev[4])
  data$X5 = mmat[5,4] * data$X4 + rnorm(size, 0, sdev[5])
  data$X6 = rnorm(size, 0, sdev[6]) 
  data$X7 = mmat[7,4] * data$X4 +  rnorm(size, 0, sdev[7])
  data$X8 = mmat[8,7] *  data$X7 + rnorm(size, 0, sdev[8]) 
  data$X9 = mmat[9,1] * data$X1 * data$X2 + mmat[9,5] *data$X5 +rnorm(size, 0, sdev[9])#interaction
  data$X10 = mmat[10,3]*data$X3 * data$X4 + mmat[10,5] * data$X5 * data$X6 + rnorm(size, 0, sdev[10])#interaction
  data$X11 = mmat[11,1] * data$X1 + mmat[11,2] * data$X2 * data$X9 + rnorm(size, 0, sdev[11]) # interaction
  data$X12= mmat[12,4] * data$X4 + mmat[12,6] * data$X6 * data$X7 + mmat[12,9] * data$X9 +
            mmat[12,10] * data$X10 + rnorm(size, 0, sdev[12]) # interaction
  
  #5 interaction terms
  return(data)
}

get.int5= function(size, mmat, sdev){
  
  data = data.frame(X1 = numeric(size))
  for(x in seq(ncol(mmat))){
    v = paste("X", x, sep = "")
    data[v] = rep(0, size)
  }  
  
  data$X1 = rnorm(size, 0, sdev[1])
  data$X2 = rnorm(size, 0, sdev[2])
  data$X3 = mmat[3,1] * data$X1 + rnorm(size, 0, sdev[3]) 
  data$X4 = rnorm(size, 0, sdev[4])
  data$X5 = 1 
  data$X6 = rnorm(size, 0, sdev[6]) 
  data$X7 = mmat[7,4] * data$X4 +  rnorm(size, 0, sdev[7])
  data$X8 = mmat[8,7] *  data$X7 + rnorm(size, 0, sdev[8]) 
  data$X9 = mmat[9,1] * data$X1 * data$X2 + mmat[9,5] *data$X5 +rnorm(size, 0, sdev[9])#interaction
  data$X10 = mmat[10,3]*data$X3 * data$X4 + mmat[10,5] * data$X5 * data$X6 + rnorm(size, 0, sdev[10])#interaction
  data$X11 = mmat[11,1] * data$X1 + mmat[11,2] * data$X2 * data$X9 + rnorm(size, 0, sdev[11]) # interaction
  data$X12= mmat[12,4] * data$X4 + mmat[12,6] * data$X6 * data$X7 + mmat[12,9] * data$X9 +
    mmat[12,10] * data$X10 + rnorm(size, 0, sdev[12]) # interaction
  
  #5 interaction terms
  return(data)
  
}
get.int6= function(size, mmat, sdev){
  
  data = data.frame(X1 = numeric(size))
  for(x in seq(ncol(mmat))){
    v = paste("X", x, sep = "")
    data[v] = rep(0, size)
  }  
  
  data$X1 = rnorm(size, 0, sdev[1])
  data$X2 = rnorm(size, 0, sdev[2])
  data$X3 = mmat[3,1] * data$X1 + rnorm(size, 0, sdev[3]) 
  data$X4 = rnorm(size, 0, sdev[4])
  data$X5 = mmat[5,4] * data$X4 + rnorm(size, 0, sdev[5])
  data$X6 = 1
  data$X7 = mmat[7,4] * data$X4 +  rnorm(size, 0, sdev[7])
  data$X8 = mmat[8,7] *  data$X7 + rnorm(size, 0, sdev[8]) 
  data$X9 = mmat[9,1] * data$X1 * data$X2 + mmat[9,5] *data$X5 +rnorm(size, 0, sdev[9])#interaction
  data$X10 = mmat[10,3]*data$X3 * data$X4 + mmat[10,5] * data$X5 * data$X6 + rnorm(size, 0, sdev[10])#interaction
  data$X11 = mmat[11,1] * data$X1 + mmat[11,2] * data$X2 * data$X9 + rnorm(size, 0, sdev[11]) # interaction
  data$X12= mmat[12,4] * data$X4 + mmat[12,6] * data$X6 * data$X7 + mmat[12,9] * data$X9 +
    mmat[12,10] * data$X10 + rnorm(size, 0, sdev[12]) # interaction
  
  #5 interaction terms
  
  return(data)
  
}

get.int10 = function(size, mmat, sdev){
  
  data = data.frame(X1 = numeric(size))
  for(x in seq(ncol(mmat))){
    v = paste("X", x, sep = "")
    data[v] = rep(0, size)
  }  
  
  data$X1 = rnorm(size, 0, sdev[1])
  data$X2 = rnorm(size, 0, sdev[2])
  data$X3 = mmat[3,1] * data$X1 + rnorm(size, 0, sdev[3]) 
  data$X4 = rnorm(size, 0, sdev[4])
  data$X5 = mmat[5,4] * data$X4 + rnorm(size, 0, sdev[5])
  data$X6 = rnorm(size, 0, sdev[6]) 
  data$X7 = mmat[7,4] * data$X4 +  rnorm(size, 0, sdev[7])
  data$X8 = mmat[8,7] *  data$X7 + rnorm(size, 0, sdev[8]) 
  data$X9 = mmat[9,1] * data$X1 * data$X2 + mmat[9,5] *data$X5 +rnorm(size, 0, sdev[9])#interaction
  data$X10 = 1
  data$X11 = mmat[11,1] * data$X1 + mmat[11,2] * data$X2 * data$X9 + rnorm(size, 0, sdev[11]) # interaction
  data$X12= mmat[12,4] * data$X4 + mmat[12,6] * data$X6 * data$X7 + mmat[12,9] * data$X9 +
    mmat[12,10] * data$X10 + rnorm(size, 0, sdev[12]) # interaction
  
  #5 interaction terms
  return(data)
  
}
get.int.data = function(mmat, sdev, size){
  
  data = get.interaction.data(size, mmat, sdev)
  
  data = rbind(data, 
               get.int3(size, mmat, sdev), 
               get.int5(size, mmat, sdev), 
               get.int6(size, mmat, sdev),
               get.int10(size, mmat, sdev))
  
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


