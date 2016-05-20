# source("sid_samplesize.R")
# #-----------------------------------------------------------------------
# # SEM:
# # X1 = N1 , N1~N(0,s1))
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
# quad = matrix(c(0,0,0,0,0,
#                 2,0,0,0,0,
#                 0,1,0,0,0,
#                 0,2,1,0,0,
#                 2,0,1,2,0),
#               nrow =5, byrow = T)
# 
# 
# size = seq(10, 2000, 2)

set.seed(6)

dat.list = lapply(size, FUN = get.pot.gaus, mmat = mmat, sdev = sdev, pot = quad)

get.int.data = function(mmat, targets, tar.values, sdev, size, pot){
  
  data = get.pot.gaus(mmat, sdev, pot, size)
  
  for( i in seq(length(targets))){
    
    int = get.int.pot(mmat, pot, targets[i], tar.values[i], sdev, size)
    data = rbind(data, int)
  }
  return(data)
}


# targets = c(2,3)
# tar.values = c(1,1)
# int.targets = list(integer(0), 2,3)
#int.size =  round(size/(length(targets) + 1))
int.data = lapply(int.size, get.int.data, targets = targets, tar.values = tar.values, mmat = mmat, sdev = sdev, pot = quad)
# int.index = lapply(size, get.int.index, targets)
# 


pdf("sid_quad_gauss.pdf")
# print("start quadratic hc")
# plot.sid(dat.list, size, hc, "hc", true.am)
# print("end hc")
# print("start quadratic gs")
# plot.sid(dat.list, size, gs, "gs", true.am)
# print("end gs")
# print("start quadratic iamb")
# plot.sid(dat.list, size, iamb, "iamb", true.am)
# print("end iamb")

print("start quadratic pc")
plot.sid(dat.list, size, pc, "pc", true.am)
print("end pc")
#plot.sid(dat.list, size, fci, "fci", mmat)
print("start quadratic ges")
plot.sid(dat.list, size, ges, "ges", true.am)
print("end ges")
print("start quadratic gies")
plot.int.sid(int.data, int.index, int.targets, method = gies, method.name = "gies", true.am, int.size)
print("end gies")

dev.off()

# qubic = matrix(c(0,0,0,0,0,
#                  3,0,0,0,0,
#                  0,1,0,0,0,
#                  0,3,1,0,0,
#                  3,0,1,3,0),
#                nrow =5, byrow = T)

dat.list = lapply(size, FUN = get.pot.gaus, mmat = mmat, sdev = sdev, pot = qubic)


# get.int.index = function(size, targets){
#   int.target.index = (rep(c(1, targets), rep(size, length(targets)+1)))
#   
# }
# 
# targets = c(2,3)
# tar.values = c(1,1)
# int.targets = list(integer(0), 2,3)
#int.size =  round(size/(length(targets) + 1))
int.data = lapply(int.size, get.int.data, targets = targets, tar.values = tar.values, mmat = mmat, sdev = sdev, pot = qubic)
# int.index = lapply(size, get.int.index, targets)

pdf("sid_qubic_gauss.pdf")
# print("start qubic hc")
# plot.sid(dat.list, size, hc, "hc", true.am)
# print("end hc")
# print("start qubic gs")
# plot.sid(dat.list, size, gs, "gs", true.am)
# print("end gs")
# print("start qubic iamb")
# plot.sid(dat.list, size, iamb, "iamb", true.am)
# print("end iamb")

print("start qubic pc")
plot.sid(dat.list, size, pc, "pc", true.am)
print("end pc")
#plot.sid(dat.list, size, fci, "fci", mmat)
print("start qubic ges")
plot.sid(dat.list, size, ges, "ges", true.am)
print("end ges")
print("start qubic gies")
plot.int.sid(int.data, int.index, int.targets, method = gies, method.name = "gies", true.am, int.size)
print("end gies")


dev.off()
