getImageDerivatives<- function(im,translation=TRUE,rotation=TRUE,skews=FALSE,scales=FALSE){
  ## spatial derivatives(x,y,z,roll,pitch,yaw, skews and scales)
  ## see http://www.johndcook.com/blog/2015/11/21/numerical-differentiation/ 
  cl<-class(im)
  if(cl!="niftiImage" && cl!="array"){
    stop(paste("class: ",cl," is not supported.\n niftiImage and arrays are supported"))
  }
  # list to hold derivatives
  A <-list()
  # small value for h
  h<-1e-3
  d<-dim(drop(im))
  if(length(d)!=3){stop("im should have 3 and only 3 dimensions")}
if(translation){  
  aff1<-buildAffine(translation=c(h,0,0),source = im)
  aff2<-buildAffine(translation=c(-h,0,0),source = im)
  xph<-applyTransform(transform = aff1,x =im)
  xmh<-applyTransform(transform = aff2,x =im)
  A[[1]]<- (xph-xmh)/(2*h)
  
  aff1<-buildAffine(translation=c(0,h,0),source = im)
  aff2<-buildAffine(translation=c(0,-h,0),source = im)
  xph<-applyTransform(transform = aff1,x = im)
  xmh<-applyTransform(transform = aff2,x = im)
  A[[2]]<- (xph-xmh)/(2*h)
  
  aff1<-buildAffine(translation=c(0,0,h),source = im)
  aff2<-buildAffine(translation=c(0,0,-h),source = im)
  xph<-applyTransform(transform = aff1,x = im)
  xmh<-applyTransform(transform = aff2,x = im)
  A[[3]]<- (xph-xmh)/(2*h)
}
if(rotation){
  aff1<-buildAffine(angles=c(h,0,0),source = im)
  aff2<-buildAffine(angles=c(-h,0,0),source = im)
  xph<-applyTransform(transform = aff1,x = im)
  xmh<-applyTransform(transform = aff2,x = im)
  A[[4]]<- (xph-xmh)/(2*h)
  
  aff1<-buildAffine(angles=c(0,h,0),source = im)
  aff2<-buildAffine(angles=c(0,-h,0),source = im)
  xph<-applyTransform(transform = aff1,x = im)
  xmh<-applyTransform(transform = aff2,x = im)
  A[[5]]<- (xph-xmh)/(2*h)
  
  aff1<-buildAffine(angles=c(0,0,h),source = im)
  aff2<-buildAffine(angles=c(0,0,-h),source = im)
  xph<-applyTransform(transform = aff1,x = im)
  xmh<-applyTransform(transform = aff2,x = im)
  A[[6]]<- (xph-xmh)/(2*h)
  }
if(skews){
    aff1<-buildAffine(skews=c(h,0,0),source = im)
    aff2<-buildAffine(skews =c(-h,0,0),source = im)
    xph<-applyTransform(transform = aff1,x = im)
    xmh<-applyTransform(transform = aff2,x = im)
    A[[7]]<- (xph-xmh)/(2*h)
    
    aff1<-buildAffine(skews=c(0,h,0),source = im)
    aff2<-buildAffine(skews=c(0,-h,0),source = im)
    xph<-applyTransform(transform = aff1,x = im)
    xmh<-applyTransform(transform = aff2,x = im)
    A[[8]]<- (xph-xmh)/(2*h)
    
    aff1<-buildAffine(skews=c(0,0,h),source = im)
    aff2<-buildAffine(skews=c(0,0,-h),source = im)
    xph<-applyTransform(transform = aff1,x = im)
    xmh<-applyTransform(transform = aff2,x = im)
    A[[9]]<- (xph-xmh)/(2*h)
  }
if(scales){
    aff1<-buildAffine(scales=c(1+h,1,1),source = im,target = im)
    aff2<-buildAffine(scales=c(1-h,1,1),source = im,target = im)
    xph<-applyTransform(transform = aff1,x = im)
    xmh<-applyTransform(transform = aff2,x = im)
    A[[10]]<- (xph-xmh)/(2*h)
    
    aff1<-buildAffine(scales=c(1,1+h,1),source = im,target = im)
    aff2<-buildAffine(scales=c(1,1-h,1),source = im,target = im)
    xph<-applyTransform(transform = aff1,x = im)
    xmh<-applyTransform(transform = aff2,x = im)
    A[[11]]<- (xph-xmh)/(2*h)
    
    aff1<-buildAffine(scales=c(1,1,1+h),source = im,target = im)
    aff2<-buildAffine(scales=c(1,1,1-h),source = im,target = im)
    xph<-applyTransform(transform = aff1,x = im)
    xmh<-applyTransform(transform = aff2,x = im)
    A[[12]]<- (xph-xmh)/(2*h)
  }
  
  return(A)
}
