# source("sid_samplesize.R")
# #-----------------------------------------------------------------------
# # SEM:
# # X1 = N1 , N1~exp(lambda1)
# # X2 = b1_2 * X1 + N2
# # X3 = b2_3 * X2 + N3
# # X4 = b2_4 * X2 + b3_4 * X3 + N4
# # X5 = b1_5 * X1 + b3_5 * X3 + b4_5 * X4 + N5
# #-----------------------------------------------------------------------
# 
# sdev = c(0.7, 0.1, 0.2, 0.8, 0.3)

# 
# mmat = matrix( c(0,    0,    0,   0,   0,
#                  0.3,  0,    0,   0,   0,
#                  0,    0.6,  0,   0,   0,
#                  0,    0.14, 0.6, 0,   0,
#                  -0.5, 0,   -0.7, 0.9, 0), 
#                nrow = 5, byrow = T)
# 
#size = seq(1990, 2000, 2)
set.seed(6)

dat.list = lapply(size, get.lin.logis.gaus, mmat = mmat, logistic = logistic, l = l, k = k, sdev = sdev)

get.int.data = function(mmat, logistic, targets, tar.values, l, k, sdev, size){
  
  data = get.lin.logis.gaus(mmat, logistic, l,k,sdev,size)
  
  for( i in seq(length(targets))){
    
    int = get.int.logis(mmat, logistic, targets[i], tar.values[i], l, k, sdev, size)
    data = rbind(data, int)
  }
  return(data)
}

int.data = lapply(int.size, get.int.data, targets = targets, tar.values = tar.values, mmat = mmat, sdev = sdev, l = l, k = k, logistic = logistic)




pdf("logistic.pdf")
# print("start logistic hc")
# plot.sid(dat.list, size, hc, "hc", true.am)
# print("end hc")
# print("start logistic gs")
# plot.sid(dat.list, size, gs, "gs", true.am)
# print("end gs")
# print("start logistic iamb")
# plot.sid(dat.list, size, iamb, "iamb", true.am)
# print("end iamb")

print("start logistic pc")
plot.sid(dat.list, size, pc, "pc", true.am)
print("end pc")
#plot.sid(dat.list, size, fci, "fci", mmat)
print("start logistic ges")
plot.sid(dat.list, size, ges, "ges", true.am)
print("end ges")
print("start logistic gies")
plot.int.sid(int.data, int.index, int.targets, method = gies, method.name = "gies", true.am, int.size)
print("end gies")

dev.off()
