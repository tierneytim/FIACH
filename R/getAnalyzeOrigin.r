getAnalyzeOrigin<-function(fname){
  
  # check if the image is actually an analyze file
  hdr<-RNifti::dumpNifti(fname)
  if(hdr$magic!=""){stop("This file does not appear to be an Analyze file")}
  
  # work out if the file is the header or the image
  hdr<-grepl(pattern = "[.]hdr",x = fname)
  img<-grepl(pattern = "[.]img",x = fname)
  
  # if the image file is passed replace the image extension with .hdr
  if(img){
    fin<-gsub(pattern ="[.]img",replacement =".hdr",x = fname)
  }else{
    fin<-fname
  }
  
  # open the file 
  fid<-gzfile(fin,"rb")
  
  # try the platform endianness to read the file
  endian<-.Platform$endian
  sizeof.hdr <- readBin(fid, integer(), size=4, endian=endian)
  
  # test of the endianess is appropriate
  if (sizeof.hdr != 348) {
  
    # if innapropriate close file and swap endian
    close(fid)
    endian <- "swap"
  
    # reopen the file   
    fid <- gzfile(fin, "rb")

    # get proper size of header and check we've got it right
    sizeof.hdr <- readBin(fid, integer(), size=4, endian=endian)
    if (sizeof.hdr != 348){stop("endianness could not be established and file could not be read")}
  
  }
  
  # read header as far as 253rd byte: note we have already read 4 bytes for the size of header
  # see http://grahamwideman.com/gw/brain/analyze/formatdoc.htm
  rubbish<- readBin(fid, raw(),n = 249 ,size=1, endian=endian)
  
  # read the 3 out of 5 bytes holding the origin
  origin<-readBin(fid, integer(), n = 3, size=2, endian=endian)
  
  #close the file
  close(fid)
  
  #return the origin
  return(origin)
} 

