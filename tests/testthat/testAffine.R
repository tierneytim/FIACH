library(Rcpp)
library(FIACH)
library(RNiftyReg)
library(mmand)
#sourceCpp("C:/Users/Nerevar/Desktop/work/applyAffine.cpp")
  Sys.setenv("PKG_CXXFLAGS"="-fopenmp")
  Sys.setenv("PKG_LIBS"="-fopenmp")
 sourceCpp("D:/applyAffine.cpp",dryRun = FALSE,verbose = TRUE,rebuild = TRUE)

n<-64
trans<-15
test <- readNifti(system.file("extdata", "flash_t1.nii.gz",package="RNiftyReg"))
pix<-pixdim(test)
d<-dim(test)


aff2<-buildAffine(translation = c(6,6,6),source = test,
                  angles = c(.05,.1,.3),
                  skews = c(.03,.4,.1),
                  scales=c(2,2,2))
pixt<-pixdim(attributes(aff2)$target)

Mt<-solve(xform(attributes(aff2)$target))
sc<-diag(4)
sc[1:3,4]<-c(-1,1,1)
Ma<-solve(aff2%*%sc)
Ms<-xform(test)
aff1<-Mt%*%Ma%*%Ms

#aff1<-sc%*%aff1


system.time(reg1<-applyAffine(test,aff1,outDim = dim(test)*2))
system.time(reg2<-applyTransform(transform = aff2,x = test,internal = NA))
dim(reg1)
dim(reg2)
l<-list()
l[[1]]<-(1:d[1])+trans/2
l[[2]]<-(1:d[2])-trans/2
l[[3]]<-(1:d[3])-trans/2
 system.time(reg3<-resample(test,points = l,kernel = mnKernel(1,0)))
# 
 viewNew(reg1-reg2[,,])
#  microbenchmark::microbenchmark(applyTransform(transform = aff2,x = test),applyAffine(test,aff1),unit = "relative")
#  microbenchmark::microbenchmark(resample(test,points = l,kernel = mnKernel(1,0)))
 

# viewR(reg1)
# viewR(reg2)
mean((reg1-reg2)^2)

viewNew(reg1[,,]-reg3)

d1<-do.call("cbind",getImageDerivatives(reg1,rotation = FALSE))
d2<-do.call("cbind",getImageDerivatives(reg2,rotation = FALSE))
dtest<-do.call("cbind",getImageDerivatives(test,rotation = FALSE))
var(dtest[!kmeansMask(test),])
var(d2[!kmeansMask(test),])
var(d1[!kmeansMask(test),])