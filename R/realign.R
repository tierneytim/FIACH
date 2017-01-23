.getOrigin<-function(fname){
  
  gzipped<-grepl(pattern = "gz$",x = fname)
  if(gzipped){stop("getOrigin does not support gzipped headers")}
  
  hdr<-grepl(pattern = "[.]hdr",x = fname)
  img<-grepl(pattern = "[.]img$",x = fname)
  if(img){
    fin<-gsub(pattern ="[.]img$",replacement =".hdr",x = fname)
  }else{
    fin<-fname
  }
  
  fid <- file(fin, "rb")
  endian<-.Platform$endian
  sizeof.hdr <- readBin(fid, integer(), size=4, endian=endian)
  if (sizeof.hdr != 348) {
    close(fid)
    endian <- "swap"
    fid <- file(fin, "rb")
  }
  rubbish<- readBin(fid, raw(),n = 249 ,size=1, endian=endian)
  origin<-readBin(fid, integer(), n = 5, size=2, endian=endian)[3:5]
  close(fid)
  return(origin)
} 

.qformFix<-function(image,origin){
  hdr<-dumpNifti(image)
  mat<-diag(4)
  d<-dim(image)[1:3]
  if(any(is.na(dim(image)[1:3]))){stop("qformFix does not support Images with less than 3 dimensions")}
  diag(mat)[1:3]<-abs(pixdim(image)[1:3])*c(-1,1,1)
  if(hdr$magic==""){
    trans<- -((origin-1)*abs(pixdim(image)[1:3]))
    mat[1:3,4]<-trans*c(-1,1,1)
    invisible(qform(image)<-structure(.Data = mat,code=2))
    invisible(sform(image)<-structure(.Data = mat,code=2))
  }
}                                                           
.makeA <- function(im){
  ## spatial derivatives(x,y,z,roll,pitch,yaw)
  ## see http://www.johndcook.com/blog/2015/11/21/numerical-differentiation/ 
  
  A <-matrix(0,nrow = length(im),ncol = 6)
  h<-1e-3
  
  aff1<-buildAffine(translation=c(h,0,0),source = im)
  aff2<-buildAffine(translation=c(-h,0,0),source = im)
  xph<-applyTransform(transform = aff1,x =im)
  xmh<-applyTransform(transform = aff2,x =im)
  A[,1]<- (xph-xmh)/(2*h)
  
  aff1<-buildAffine(translation=c(0,h,0),source = im)
  aff2<-buildAffine(translation=c(0,-h,0),source = im)
  xph<-applyTransform(transform = aff1,x = im)
  xmh<-applyTransform(transform = aff2,x = im)
  A[,2]<- (xph-xmh)/(2*h)
  
  aff1<-buildAffine(translation=c(0,0,h),source = im)
  aff2<-buildAffine(translation=c(0,0,-h),source = im)
  xph<-applyTransform(transform = aff1,x = im)
  xmh<-applyTransform(transform = aff2,x = im)
  A[,3]<- (xph-xmh)/(2*h)
  
  aff1<-buildAffine(angles=c(h,0,0),source = im)
  aff2<-buildAffine(angles=c(-h,0,0),source = im)
  xph<-applyTransform(transform = aff1,x = im)
  xmh<-applyTransform(transform = aff2,x = im)
  A[,4]<- (xph-xmh)/(2*h)
  
  aff1<-buildAffine(angles=c(0,h,0),source = im)
  aff2<-buildAffine(angles=c(0,-h,0),source = im)
  xph<-applyTransform(transform = aff1,x = im)
  xmh<-applyTransform(transform = aff2,x = im)
  A[,5]<- (xph-xmh)/(2*h)
  
  aff1<-buildAffine(angles=c(0,0,h),source = im)
  aff2<-buildAffine(angles=c(0,0,-h),source = im)
  xph<-applyTransform(transform = aff1,x = im)
  xmh<-applyTransform(transform = aff2,x = im)
  A[,6]<- (xph-xmh)/(2*h)
  
  return(A)
}
.inVox <- function(A0,b,quality=.9,plotConvergence=FALSE){
  rat<-1
  Alpha<-cbind(A0,b)
  cps<-apply(Alpha,MARGIN = 1,tcrossprod)
  Alpha<-crossprod(Alpha)
  xyz<-which(is.finite(b),arr.ind = TRUE)
  det0<-det(Alpha)
  det1<-det0
  det1/det0
  
  
  while(det1/det0>=quality){
    dets<-.allDets(cps = cps,Alpha = Alpha)
    msk<-order(det1-dets)
    msk<-msk[1:round(length(dets)/10)]
    A0<-A0[-msk,]
    b<-b[-msk]
    cps<-cps[,-msk]
    xyz<-xyz[-msk,]
    Alpha<-cbind(A0,b)
    Alpha<-crossprod(Alpha)
    det1<-det(Alpha)
    rat<-c(rat,det1/det0)
    if(plotConvergence){ts.plot(rat,ylab = "Quality",xlab = "Iteration",main="Identifying Interesting Voxels")}
  }
  return(list(xyz,A0,b))
}        
realign<-function(files,quality=.9,write=TRUE,plotConvergence=FALSE){
  ##########################
  ##### ARG checking #######
  ##########################
  print("Argument Check")
  if(!is.character(files)){stop("input should be character strings")}
  exists<-file.exists(files)
  if(!all(exists)){stop("At least one of the specified functional files does not exist")}
  if(!is.numeric(quality)){stop("quality should be a real number between 0 and 1")}
  if(quality==0){stop("quality must be a greater than 0")}
  if(quality>1){stop("quality must be less than or equal to 1")}
  dims<-do.call("rbind",lapply(files,function(x){dumpNifti(x)$dim}))
  if(any(dims[,5]>1)){stop("Currently this function does not support  4D files")}
  if(any(dims[1,]!=t(dims))){stop("Not all image dimensions are the same as the target image")}
  if(!is.logical(write)){stop("write must be TRUE or FALSE")}
  if(!is.logical(plotConvergence)){stop("plotConvergence must be TRUE or FALSE")}
  ##########################
  ### FILE READ ############
  ##########################
  print("Assuming target image is first image")
  targ<-RNiftyReg::readNifti(files[1])
  anlz<-dumpNifti(files[1])$magic==""
  qf<-xform(targ)
  pix<-pixdim(targ)[1:3]
  origin<-(abs(qf[1:3,4])+abs(pix))/pix
  if(anlz){origin<-.getOrigin(files[1])}
  .qformFix(targ,origin)
  print("Fixing target q/s form")
  pix <- pixdim(targ)                                                                 ## pixel dimensions(assume all are same)
  aff0<-buildAffine(scales = pix[1:3]/4,source = drop(targ))                          ## initial affine used to rescale  target 
  print("Resampling target to Lower Resolution")
  b <- applyTransform(drop(targ),transform = aff0,internal=NA)                                    ## (Ax = b)
  ##########################
  ###### FILE MANAGEMENT ###
  ##########################
  print("Creating output filenames")
  dir<-dirname(files)
  fname<-basename(files)
  outname<-paste(dir,"/r",fname,sep="")
  rpOutName<-paste(dir[1],"rp.txt",sep="/")
  rpPlotOutName<-paste(dir[1],"rp.pdf",sep="/")
  outCode<-getDatatype(input = files)
  unsigned<-(outCode=="char"||outCode=="uint16"||outCode=="uint32"||outCode=="uint64")
  if(unsigned){print("Output will have  negative values set to zero as input has an unsigned datatype")}
  ##########################
  ####### GRADIENTS ########
  ##########################
  print("Constructing Derivatives")
  A0 <- .makeA(b)                                                                           ## basis set of derivatives
  ###########################
  ### INTERESTING VOXELS ####
  ###########################
  print("Identifying Interesting Voxels")
  co <- .inVox(A0,b,quality = quality,plotConvergence = plotConvergence)                                                                        ## interesting voxels
  ###########################
  ###### REG VARIABLES ######
  ###########################
  print("Initialising Variables for Registration")
  qA <- qr(co[[2]])                                                                        ## decompose the masked basis set for linear least squares using QR
  ss<-Inf                                                                                   ## initialize sum of squares
  maxit<-64                                                                                ## max iterations
  reg<-b
  rp<-matrix(0,nrow=sum(dims[,5]),ncol = 6)                                                ## declare matrix to hold realignment paramters
  aff<-buildAffine(source = b)                                                             ## declare affine to be used in registration
  bsub <- co[[3]]                                                                          ## masked target as a vector
  sbsub<-sum(bsub)                                                                         ## sum of masked target as a vector(prevent recalculation in main loop)
  ############################
  ##### REGISTRATION #########
  ############################
  print("Registering and Reslicing")

  pb <- txtProgressBar(style = 3)
  for (j in 1:length(files)) {                                                               ## main loop over images
    if(j==1){
      reg<-targ
      if(write){writeNifti(image = reg,file = outname[j],datatype = outCode)}
      rp[j,]<-0
      setTxtProgressBar(pb, j/length(files))
      next
    }
    sourceF<-readNii(input = files[j])
    .qformFix(sourceF,origin = origin)
    rfunc<-applyTransform(transform = aff0,interpolation = 1,x = sourceF,internal=NA)
    for (i in 1:maxit) {                                                                    ## registration loop
      Fl <- applyTransform(transform = aff,x = rfunc,interpolation = 1,internal=NA)[co[[1]]]              ## apply potentially updated transform                                                       ## the masked target as a vector
      sc <- sbsub/sum(Fl)                                                                     ## ratio of the source and the target 
      b1 <- bsub - Fl * sc                                                                    ## scale source so mean intensity ~= target
      pss<-ss;                                                                                ## declare the previous sum of squares
      ss<-sum(b1*b1)                                                                          ## calculate new sum of squares
      if(i==1){conv<-Inf}else{conv<-(pss-ss)/pss}                                             ## define the convergence criteria
      if(pss == 0){break()}
      if ( conv< 1e-8 || i == maxit){break()}                                                 ## finish if convergerd or max iterations are reached
      soln<- qr.coef(qr = qA,y = b1)                                                          ## use base qr to get coef 
      tmp<-buildAffine(translation = soln[1:3],angles = soln[4:6],source = rfunc)             ## build affine with updates
      aff<-composeTransforms(transform1 = tmp,aff)                                            ## compose updates with previous solution
    }                                                                                       ## end registraion loop          
    hmm<-decomposeAffine(aff)                                                                 ## decompose the rigid transformation
    rp[j,]<-c(hmm$translation,hmm$angles)                                                     ## get translations and rotations
    finAff<-buildAffine(translation = rp[j,1:3],angles = rp[j,4:6],source = sourceF)          ## build the final affine
    reg<-applyTransform(transform = finAff,x = sourceF,internal = NA)                                       ## apply final affine to input images
    if(unsigned){
     reg[reg<0]<-0
      }
    if(write){writeNifti(image = reg,file = outname[j],datatype = outCode)}
    setTxtProgressBar(pb, j/length(files))
  } 
  close(pb)
  #############################
  ######## PARAMETRS ##########
  #############################
  if(write){
    pdf(file = rpPlotOutName,width = 8.27,height = 11.69)
    par(mfrow=c(2,1))
    ts.plot(rp[,1:3,drop=FALSE],main="Translation Paramteres",xlab="Time (scans)",ylab="Distance (mm)")
    ts.plot(rp[,4:6,drop=FALSE]*57,main="Rotation Parameters",xlab="Time (scans)",ylab="Rotation (Degrees)")
    dev.off()
    write.table(rp,file = rpOutName,row.names = FALSE,col.names = FALSE)
  }
  return(rp)
}
