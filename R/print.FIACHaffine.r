print.FIACHAffine <- function(x, ...) {
  affineClass <- class(x)
  if (affineClass != "FIACHAffine") {
    stop(paste("class: ", affineClass, " is not supported.\n objects of class FIACHAffine are supported "))
}  
  cat("FIACH affine matrix:\n")
  lines <- apply(format(x,scientific=FALSE), 1, paste, collapse="  ")
  cat(paste(lines, collapse="\n"))
}
