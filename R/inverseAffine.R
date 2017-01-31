inverseAffine<-function(aff){
  
  affineClass<-class(aff)
  
  if(affineClass!="FIACH Affine"){
    stop(paste("class: ",affineClass," is not supported.\n objects of class FIACH Affine are supported "))
  }
  sx<-attr(aff,"sourceXform")  
  tx<-attr(aff,"targetXform")
  tDim<-attr(aff,"targetDim")
  sDim<-attr(aff,"sourceDim")
  
  out<-solve(aff)
  attr(out,"targetXform")<-sx
  attr(out,"sourceXform")<-tx
  attr(out,"targetDim")<-sDim
  attr(out,"sourceDim")<-tDim
  class(out)<-"FIACH Affine"
  
  return(out)
}
