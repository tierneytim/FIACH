print.FIACHAffine <- function(x, ...) {
  cat("FIACH affine matrix:\n")
  lines <- apply(format(x,scientific=FALSE), 1, paste, collapse="  ")
  cat(paste(lines, collapse="\n"))
}
