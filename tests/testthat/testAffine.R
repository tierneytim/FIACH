library(FIACH)
library(RNiftyReg)
################################
###### source---> target #######
################################
source <- readNifti(system.file("extdata", "epi_t2.nii.gz",package="RNiftyReg"))
target <- readNifti(system.file("extdata", "flash_t1.nii.gz",package="RNiftyReg"))

affReg<-buildAffine(source = source,target = target)
affTim<-createAffine(source=source,target=target)

reg1<-applyTransform(affReg,x = source)
reg2<-applyAffine(source=source,aff = affTim)
sd(reg1-reg2)

invAffReg<-invertAffine(affReg)
invAffTim<-inverseAffine(affTim)

reg3<-applyTransform(invAffReg,x = reg1)                  
var(reg3-source)
reg4<-applyAffine(source = reg2,aff = invAffTim)
var(reg4-source)
###############################
# source---> source(rigid) ####
###############################
source<-target

affReg<-buildAffine(source = source,
                    translation = c(6.4,-3.5,9.4),
                    angles=c(.1,-.2,.3)
                    )
affTim<-createAffine(source=source,
                     translation = c(6.4,-3.5,9.4),
                     angles=c(.1,-.2,.3)
                     )

reg1<-applyTransform(affReg,x = source)
reg2<-applyAffine(source=source,aff = affTim)
sd(reg1-reg2)

invAffReg<-invertAffine(affReg)
invAffTim<-inverseAffine(affTim)

reg3<-applyTransform(invAffReg,x = reg1)                  
var(reg3-source)
reg4<-applyAffineTim(source = reg2,aff = invAffTim)
var(reg4-source)
###############################
# source---> source(affine) ###
###############################
translation =c(0,0,0)
angles = c(0,0,0)
skews = c(0,0,0)
scales = c(1.2,1.3,1.7)


affReg<-buildAffine(source = source,
                      translation =translation,
                      angles = angles,
                      skews = skews,
                      scales = scales)
  
affTim<-createAffine(source = source,
                    translation =translation,
                    angles = angles,
                    skews = skews,
                    scales = scales)

reg1<-applyTransform(affReg,x = source)
reg2<-applyAffine(source=source,aff = affTim)
sd(reg1-reg2)

invAffReg<-invertAffine(affReg)
invAffTim<-inverseAffine(affTim)

reg3<-applyTransform(invAffReg,x = reg1)                  
var(reg3-source)
reg4<-applyAffine(source = reg2,aff = invAffTim)
var(reg4-source)
###############################
# source---> target(rigid)  ###
###############################
source <- readNifti(system.file("extdata", "epi_t2.nii.gz",package="RNiftyReg"))
target <- readNifti(system.file("extdata", "flash_t1.nii.gz",package="RNiftyReg"))


translation =c(6.4,-3.1,9.34)
angles = c(-.01,.3,.1)
skews = c(0,0,0)
scales = c(1,1,1)


affReg<-buildAffine(source = source,
                    translation =translation,
                    angles = angles,
                    skews = skews,
                    scales = scales,target = target)

affTim<-createAffine(source = source,
                     translation =translation,
                     angles = angles,
                     skews = skews,
                     scales = scales,target=target)

reg1<-applyTransform(affReg,x = source)
reg2<-applyAffine(source=source,aff = affTim)
sd(reg1-reg2)

invAffReg<-invertAffine(affReg)
invAffTim<-inverseAffine(affTim)

reg3<-applyTransform(invAffReg,x = reg1)                  
var(reg3-source)
reg4<-applyAffine(source = reg2,aff = invAffTim)
var(reg4-source)
###############################
# source---> target(affine) ###
###############################
source <- readNifti(system.file("extdata", "epi_t2.nii.gz",package="RNiftyReg"))
target <- readNifti(system.file("extdata", "flash_t1.nii.gz",package="RNiftyReg"))


translation =c(6.4,-3.1,9.34)
angles = c(-.01,.3,.1)
skews = c(.1,.03,-.2)
scales = c(1.13,1.12,1.17)


affReg<-buildAffine(source = source,
                    translation =translation,
                    angles = angles,
                    skews = skews,
                    scales = scales,target = target)

affTim<-createAffine(source = source,
                     translation =translation,
                     angles = angles,
                     skews = skews,
                     scales = scales,target=target)

reg1<-applyTransform(affReg,x = source)
reg2<-applyAffine(source=source,aff = affTim)
sd(reg1-reg2)

invAffReg<-invertAffine(affReg)
invAffTim<-inverseAffine(affTim)

reg3<-applyTransform(invAffReg,x = reg1)                  
var(reg3-source)
reg4<-applyAffine(source = reg2,aff = invAffTim)
var(reg4-source)



###############################
##### DECOMPOSITION ###########
###############################



invAff<-affTim
trans<-invAff[1:3,4]
subAff<-invAff[1:3,1:3]
decomp<-chol(crossprod(subAff))
scales<-diag(diag(decomp))
skew<-decomp%*%solve(scales)
rotation<-subAff%*%solve(decomp)

