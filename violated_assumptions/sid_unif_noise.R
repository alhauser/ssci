# source("sid_samplesize.R")
# #-----------------------------------------------------------------------
# # SEM:
# # X1 = N1 , N1~unif(a1,b1)
# # X2 = b1_2 * X1 + N2
# # X3 = b2_3 * X2 + N3
# # X4 = b2_4 * X2 + b3_4 * X3 + N4
# # X5 = b1_5 * X1 + b3_5 * X3 + b4_5 * X4 + N5
# #-----------------------------------------------------------------------
# 
# sdev = c(0.7, 0.1, 0.2, 0.8, 0.3)
a = -sqrt(3) * sdev
b = sqrt(3) * sdev
range = data.frame(a,b)# this results error terms with the same variatoin as with a normal distr and mean 0
                       # since var(u) = 1/12(b-a)

# mmat = matrix( c(0,    0,    0,   0,   0,
#                  0.3,  0,    0,   0,   0,
#                  0,    0.6,  0,   0,   0,
#                  0,    0.14, 0.6, 0,   0,
#                  -0.5, 0,   -0.7, 0.9, 0), 
#                nrow = 5, byrow = T)
# 
# size = seq(10, 2000, 2)
set.seed(6)

dat.list = lapply(size, FUN = get.lin.unif, mmat = mmat, rnge =  range)

get.int.data = function(mmat, targets, tar.values, sdev, size){
  
  data = get.lin.unif(mmat, range, size)
  
  for( i in seq(length(targets))){
    
    int = get.int.unif(mmat, targets[i], tar.values[i], range, size)
    data = rbind(data, int)
  }
  return(data)
}

# get.int.index = function(size, targets){
#   int.target.index = (rep(c(1, targets), rep(size, length(targets)+1)))
#   
# }
# 
# targets = c(2,3)
# tar.values = c(1,1)
# int.targets = list(integer(0), 2,3)
int.data = lapply(int.size, get.int.data, targets = targets, tar.values = tar.values, mmat = mmat, sdev = sdev)
#int.index = lapply(size, get.int.index, targets)



pdf("sid_unif_noise.pdf")
# print("start unif noise hc")
# plot.sid(dat.list, size, hc, "hc", true.am)
# print("end hc")
# print("start unif noise gs")
# plot.sid(dat.list, size, gs, "gs", true.am)
# print("end gs")
# print("start unif noise iamb")
# plot.sid(dat.list, size, iamb, "iamb", true.am)
# print("end iamb")

print("start unif noise pc")
plot.sid(dat.list, size, pc, "pc", true.am)
print("end pc")
#plot.sid(dat.list, size, fci, "fci", mmat)
print("start unif noise ges")
plot.sid(dat.list, size, ges, "ges", true.am)
print("end ges")
print("start unif noise gies")
plot.int.sid(int.data, int.index, int.targets, method = gies, method.name = "gies", true.am, int.size)
print("end gies")
dev.off()
