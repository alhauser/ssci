########################################################################
# this scribt is called to by a scribt named semX.R  in a subfolder 
# which defines a sem.
# thsi scribt then calles the lin gauss data and calles the function in
# sid_samplesize.R to create the sid plots
########################################################################


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

set.seed(6)

dat.list = lapply(size, FUN = get.lin.gaus, mmat = mmat, sdev = sdev, sem = NULL)

get.int.data = function(mmat, targets, tar.values, sdev, size){
  
  data = get.lin.gaus(mmat = mmat, sdev = sdev, size = size)
  
  for( i in seq(length(targets))){
    
    int = get.int.lin.gaus(mmat, targets[i], tar.values[i], sdev, size)
    data = rbind(data, int)
  }
  return(data)
}


# targets = c(2,3)
# tar.values = c(1,1)
# int.targets = list(integer(0), 2,3)
#int.size =  round(size/(length(targets) + 1))
int.data = lapply(int.size, get.int.data, targets = targets, tar.values = tar.values, mmat = mmat, sdev = sdev)


pdf("sid_lin_gauss.pdf")
# print("start lin gauss hc")
# plot.sid(dat.list, size, hc, "hc", true.am)
# print("end hc")
# print("start lin gauss gs")
# plot.sid(dat.list, size, gs, "gs", true.am)
# print("end gs")
# print("start lin gauss iamb")
# plot.sid(dat.list, size, iamb, "iamb", true.am)
# print("end iamb")

print("start lin gauss pc")
plot.sid(dat.list, size, pc, "pc", true.am)
print("end pc")
#plot.sid(dat.list, size, fci, "fci", mmat)
print("start lin gauss ges")
plot.sid(dat.list, size, ges, "ges", true.am)
print("end ges")
print("start lin gauss gies")
plot.int.sid(int.data, int.index, int.targets, method = gies, method.name = "gies", true.am, int.size)
print("end gies")


dev.off()
