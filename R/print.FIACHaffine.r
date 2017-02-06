print.FIACHAffine <- function(x, ...) {
  cat("FIACH affine matrix:\n")
  lines <- apply(format(x,scientific=FALSE), 1, paste, collapse="  ")
  cat(paste(lines, collapse="\n"))
}


file<-system.file("extdata","mni.nii.gz",package="FIACH")
source <- readNii(file)

translation =c(6.4,-3.1,9.34)
angles = c(-.01,.3,.1)
skews = c(.1,.03,-.2)
scales = c(1.13,1.12,1.17)


affReg<-createAffine(source = source,
                     translation =translation,
                     angles = angles,
                     skews = skews,
                     scales = scales)