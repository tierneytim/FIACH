fiachMerge<-function(input,compress=TRUE){
  outCode<-getDatatype(input = input)
  data<-readNii(input = input,fourD=TRUE)
  
  dir<-dirname(input)[1]
  base<-basename(input)
  basNoExt<-gsub("([.]nii|[.]img|[.]hdr|[.]nii[.]gz)",replacement = "",x = base)[1]
  if(compress){
    outname<-paste(dir,"/","4D_",basNoExt,".nii.gz",sep = "")
  }else{
    outname<-paste(dir,"/","4D_",basNoExt,sep = "")
  }
  
  RNiftyReg::writeNifti(image = data,file = outname,datatype = outCode)
}