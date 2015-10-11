reslice<-function(data,newDim,currentDim=NULL){
  dOld<-dim(data)
  if(length(dOld)==4){targ<-updateNifti(data[,,,1],template = data);dOld<-dOld[1:3]}else{
    targ<-data
  }
  if(length(newDim)!=length(dOld)){stop("No of dimensions in newDim and data differ")}
  if(!is.null(currentDim)&&length(currentDim)!=length(dOld)){stop("If currentDim is not NULL length(currnetDim) must equal length(dim(data))")}
  if(!isImage(object = data,unsure = FALSE)&& is.null(currentDim)){
    stop("Please supply current dimensions")
  }
  if(isImage(object = data,unsure = FALSE)){currentDim<-RNiftyReg::pixdim(data)[1:3]}
  
  aff<-RNiftyReg::buildAffine(scales=currentDim/newDim,source=targ)
  res<-RNiftyReg::applyTransform(transform = aff,x = data)
  if(!isImage(object = data,unsure = FALSE)){pixdim(res)<-newDim}
  return(res)
}
