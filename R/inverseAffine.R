inverseAffine<-function(aff){
  
  affineClass<-class(aff)
  
  if(affineClass!="FIACHAffine"){
    stop(paste("class: ",affineClass," is not supported.\n objects of class FIACHAffine are supported "))
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
  attr(out,"targetPixdim")<-attr(aff,"sourcePixdim")
  attr(out,"sourcePixdim")<-attr(aff,"targetixdim")
  class(out)<-"FIACHAffine"
  
  return(out)
}
