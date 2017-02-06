fiachSplit<-function(input){
  outCode<-getDatatype(input = input)
  data<-readNii(input,fourD=FALSE)
  t<-length(data)
  dir<-dirname(input)
  base<-basename(input)
  
  basNoExt<-gsub("([.]nii|[.]img|[.]hdr|[.]nii[.]gz)",replacement = "",x = base)
  num<-formatC(1:t, width = nchar(as.character(t)), format = "d", flag = "0")
  outnames<-paste(dir,"/",basNoExt,"_",num,sep = "")
  for(i in 1:t){
    RNifti::writeNifti(image = data[[i]],datatype = outCode,file = outnames[i])
  }
}