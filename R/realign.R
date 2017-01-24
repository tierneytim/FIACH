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
  
  ## files check (character and existance)
  if(!is.character(files)){stop("input should be character strings")}
  exists<-file.exists(files)
  if(!all(exists)){stop("At least one of the specified functional files does not exist")}
  
  # quality checks (numeric and between 0 and 1)
  if(!is.numeric(quality)){stop("quality should be a real number between 0 and 1")}
  if(quality==0){stop("quality must be greater than 0")}
  if(quality>1){stop("quality must be less than or equal to 1")}
  
  # dimension checks( equality and 4d files currently not allowed: i hope to fix this if i can access sub-bricks in RNifti )
  dims<-do.call("rbind",lapply(files,function(x){dumpNifti(x)$dim}))
  if(any(dims[,5]>1)){stop("Currently this function does not support  4D files")}
  if(any(dims[1,]!=t(dims))){stop("Not all image dimensions are the same as the target image")}
  
  # write check (must be logical)
  if(!is.logical(write)){stop("write must be TRUE or FALSE")}
  
  # plot convergence check (must be logical)
  if(!is.logical(plotConvergence)){stop("plotConvergence must be TRUE or FALSE")}
  ##########################
  ### FILE READ ############
  ##########################
  print("Assuming target image is first image")
  
  targ<-RNiftyReg::readNifti(files[1])
  
  # check for analyze format
  anlz<-dumpNifti(files[1])$magic==""
  
  #get the pixel dimensions
  pix<-pixdim(targ)[1:3]
  
  # if the image is analyze then try and find and origin and set an appropriate qform and signal warning.
  if(anlz){
    warning("This file appears to be an Analyze file. Orientation information may not be correct.")
    origin<-getAnalyzeOrigin(files[1])
    .qformFix(targ,origin)
    print("Fixing target q/s form")
  }
  ## affine matrix to downsample image to 4x4x4mm (same as SPM default).
  aff0<-buildAffine(scales = pix[1:3]/4,source = drop(targ))                          
  print("Resampling target to Lower Resolution")
  b <- applyTransform(drop(targ),transform = aff0,internal=NA)                                   
  ##########################
  ###### FILE MANAGEMENT ###
  ##########################
  print("Creating output filenames")
  
  # file paths into directories and files
  dir<-dirname(files)
  fname<-basename(files)
  
  # prepend filename with r
  outname<-paste(dir,"/r",fname,sep="")
  rpOutName<-paste(dir[1],"rp.txt",sep="/")
  rpPlotOutName<-paste(dir[1],"rp.pdf",sep="/")
  
  # find datatype and check whether it is unsigned
  outCode<-getDatatype(input = files)
  unsigned<-(outCode=="char"||outCode=="uint16"||outCode=="uint32"||outCode=="uint64")
  if(unsigned){print("Output will have  negative values set to zero as input has an unsigned datatype")}
  
  ##########################
  ####### GRADIENTS ########
  ##########################
  print("Constructing Derivatives")
  
  #convert the  list of 6 images to a pixel x 6 matrix
  A0 <- do.call("cbind",getImageDerivatives(b))                                                                          
  ###########################
  ### INTERESTING VOXELS ####
  ###########################
  print("Identifying Interesting Voxels")
  
  # remove voxels that contribute little to the spatial covariance matrix
  # returns masked coordinates, derivatives and target image 
  co <- .inVox(A0,b,quality = quality,plotConvergence = plotConvergence)                                                                        
  ###########################
  ###### REG VARIABLES ######
  ###########################
  print("Initialising Variables for Registration")
  
  ## decompose the masked derivatives for linear least squares using QR decomposition
  ## this is internally the default method for lm and lm.fit which compromises speed for accuracy.
  qA <- qr(co[[2]])                                                                        
  
  ## initialize sum of squares
  ss<-Inf                                                                                   
  
  ## max iterations: set to SPM default
  maxit<-64
  
  ## declare matrix to hold realignment paramters
  rp<-matrix(0,nrow=sum(dims[,5]),ncol = 6)
  
  ## masked target image  as a vector
  bsub <- co[[3]]                                                                          
  
  ## sum of masked target as a vector used to scale images in case of gross global intensity difference
  sbsub<-sum(bsub)
  ############################
  ##### REGISTRATION #########
  ############################
  print("Registering and Reslicing")

  # initialise a progress bar
  pb <- txtProgressBar(style = 3)

  # loop over j files
  for (j in 1:length(files)) {                                                               
    
    # if its the first image 
    if(j==1){
     # the target image is left unchanged
      reg<-targ
      # and written to disk
      if(write){writeNifti(image = reg,file = outname[j],datatype = outCode)}
      # the realignmet parameter matrix initial row is set to 0 because no registration occured
      rp[j,]<-0
      # update the progress bar
      setTxtProgressBar(pb, j/length(files))
      next
    }
    
    # read in the jth source  image 
    sourceF<-readNii(input = files[j])
    
    # if analyze image is detected
    if(anlz){
    # set the qform to something hopefully more appropriate
    # a warning should have triggred earlier so I won't repeat
      .qformFix(sourceF,origin = origin)
    }
    
    ## initial affine transform used to rescale  source 
    aff0<-buildAffine(scales = pix[1:3]/4,source = drop(sourceF))                          
    rfunc<-applyTransform(transform = aff0,interpolation = 1,x = sourceF,internal=NA)
    
    # identity affine that will be recursively updated sending source image into target space
    aff<-buildAffine(source = rfunc)
    
    ## registration loop
    for (i in 1:maxit) {
      # sample the image at the coordinates held in co[[1]] having applied the potentiall updated aff Transform
      Fl <- applyTransform(transform = aff,x = rfunc,interpolation = 1,internal=NA)[co[[1]]]
      
      ## ratio of the source and the target to be used as a scale factor to account for gross global differences in intensity
      sc <- sbsub/sum(Fl)                                                                     
      
      ##  difference between tareget and  scaled source image(~= mean intensity due to scaling)
      b1 <- bsub - Fl * sc                                                                    
      
      ## declare the previous sum of squares
      pss<-ss 
      
      ## calculate new sum of squares
      ss<-sum(b1*b1)                                                                          
      
      ## define the convergence criteria
      if(i==1){conv<-Inf}else{conv<-(pss-ss)/pss}
      if(pss == 0){break()}
    
      ## finish if convergerd or max iterations are reached
      if ( conv< 1e-8 || i == maxit){break()}                                                 
      
      ## the realignment parameters are obtained by solving Ax=b(qAx=b1)
      soln<- qr.coef(qr = qA,y = b1)                                                          
      
      ## construct temporary affine that contains the solution
      tmp<-buildAffine(translation = soln[1:3],angles = soln[4:6],source = rfunc)             
      
      # update final affine with the temporary one by composing them
      aff<-composeTransforms(transform1 = tmp,aff)                                            
    }                                                                                       
    ## decompose the rigid transformation and place in the realignmet parameter matrix
    hmm<-decomposeAffine(aff)                                                                 
    rp[j,]<-c(hmm$translation,hmm$angles)                                                     
    
    ## build the final affine and apply the transform
    finAff<-buildAffine(translation = rp[j,1:3],angles = rp[j,4:6],source = sourceF)          
    reg<-applyTransform(transform = finAff,x = sourceF,internal = NA)                 
    
    ## set negative values to 0 if unsigned data is being used(really should change this behaviour)
    if(unsigned){
     reg[reg<0]<-0
      }
    
    ## If writing is requested then write the resliced images to a file
    if(write){writeNifti(image = reg,file = outname[j],datatype = outCode)}
    
    # update  the progress bar
    setTxtProgressBar(pb, j/length(files))
  } 
  # close the progress bar
  close(pb)
  #############################
  ######## PARAMETRS ##########
  #############################
  # if writing is requested
  if(write){
    
    ## create a pdf of the realignemnt paramters
    pdf(file = rpPlotOutName,width = 8.27,height = 11.69)
    par(mfrow=c(2,1))
    ts.plot(rp[,1:3,drop=FALSE],main="Translation Paramteres",xlab="Time (scans)",ylab="Distance (mm)")
    ts.plot(rp[,4:6,drop=FALSE]*57,main="Rotation Parameters",xlab="Time (scans)",ylab="Rotation (Degrees)")
    dev.off()
    
    # and write the paramters to a .txt file.
    write.table(rp,file = rpOutName,row.names = FALSE,col.names = FALSE)
  }
  
  # return the parameters
  return(rp)
}
