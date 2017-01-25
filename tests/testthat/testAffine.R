library(Rcpp)
library(FIACH)
library(RNiftyReg)
library(mmand)
sourceCpp("C:/Users/Nerevar/Desktop/work/applyAffine.cpp")

n<-64
trans<-6.4
test <- readNifti(system.file("extdata", "flash_t1.nii.gz",package="RNiftyReg"))
pix<-pixdim(test)
d<-dim(test)

aff1<-buildAffine(translation = c(trans,-trans,-trans)/pix,source = test)
aff2<-buildAffine(translation = c(trans,trans,trans),source = test)

system.time(reg1<-applyAffine(test,aff1))
system.time(reg2<-applyTransform(transform = aff2,x = test))
l<-list()
l[[1]]<-(1:d[1])+trans/2
l[[2]]<-(1:d[2])-trans/2
l[[3]]<-(1:d[3])-trans/2
system.time(reg3<-resample(test,points = l,kernel = mnKernel(1,0)))

viewNew(reg1-reg2)
# microbenchmark::microbenchmark(applyTransform(transform = aff2,x = test))
# microbenchmark::microbenchmark(applyAffine(test,aff1))


# viewR(reg1)
# viewR(reg2)

viewNew(reg1-reg3)


