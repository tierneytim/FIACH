createAffine<-function (translation = c(0,0,0), scales = c(1,1,1), skews = c(0,0,0), angles = c(0,0,0), source = NULL, target = NULL){
  
  if(length(translation)!=3 || !is.numeric(translation)){
    stop("translation must be a numeric vector of length 3")
  }
  
  if(length(scales)!=3 || !is.numeric(scales)){
    stop("scales must be a numeric vector of length 3")
  }
  
  if(length(skews)!=3 || !is.numeric(skews)){
    stop("skews must be a numeric vector of length 3")
  }
  
  if(length(angles)!=3 || !is.numeric(angles)){
    stop("angles must be a numeric vector of length 3")
  }
  
  sourceClass<-class(source)
  targetClass<-class(source)
  
  if(sourceClass!="niftiImage"  && sourceClass!="array"){
    stop(paste("class: ",sourceClass," is not supported.\n niftiImage and arrays are supported "))
  }
  if(targetClass!="niftiImage"  && targetClass!="array" && targetClass!="NULL"){
    stop(paste("class: ",targetClass," is not supported.\n niftiImage,arrays and objects of class NULL are supported "))
  }
  if(any(scales<=0)){
    stop("all scales hould be greater than zero")
  }
  
  
  
  sourceXform<-xform(source)
  sourceDim<-dim(source)
  
  if(is.null(target)){
    targetXform<-sourceXform
    targetDim<-floor(scales*dim(source))
    scaleMat<-diag(scales)
    targetXform[1:3,1:3]<-solve(scaleMat)%*%targetXform[1:3,1:3]
    scales<-c(1,1,1)
  }else{
    targetXform<-xform(target)
    targetDim<-dim(target)
  }
  affine <- diag(4)
  
  rotationX <- rotationY <- rotationZ <- skewMatrix <- diag(3)
  cosAngles <- cos(angles)
  sinAngles <- sin(angles)
  rotationX[2:3,2:3] <- c(cosAngles[1], -sinAngles[1], sinAngles[1], cosAngles[1])
  rotationY[c(1,3),c(1,3)] <- c(cosAngles[2], sinAngles[2], -sinAngles[2], cosAngles[2])
  rotationZ[1:2,1:2] <- c(cosAngles[3], -sinAngles[3], sinAngles[3], cosAngles[3])
  skewMatrix[c(4,7,8)] <- skews
  
  affine[1:3,1:3] <- rotationX %*% rotationY %*% rotationZ %*% skewMatrix %*% diag(scales)
  affine[1:3,4] <- translation
  attr(affine,"sourceXform")<-sourceXform
  attr(affine,"targetXform")<-targetXform
  attr(affine,"sourceDim")<-sourceDim
  attr(affine,"targetDim")<-targetDim
  class(affine)<-"FIACH Affine"
  return(affine)
} 
