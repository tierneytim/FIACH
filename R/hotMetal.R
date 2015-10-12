hotMetal<-function(n){
  orig <- c("#010000", "#0C0000", "#170000", "#210000", "#2C0000", 
            "#360000", "#410000", "#4C0000", "#560000", "#610000", 
            "#6C0000", "#760000", "#810000", "#8B0000", "#960000", 
            "#A10000", "#AB0000", "#B60000", "#C10000", "#CB0000", 
            "#D60000", "#E00000", "#EB0000", "#F60000", "#FF0100", 
            "#FF0C00", "#FF1700", "#FF2100", "#FF2C00", "#FF3600", 
            "#FF4100", "#FF4C00", "#FF5600", "#FF6100", "#FF6C00", 
            "#FF7600", "#FF8100", "#FF8B00", "#FF9600", "#FFA100", 
            "#FFAB00", "#FFB600", "#FFC100", "#FFCB00", "#FFD600", 
            "#FFE000", "#FFEB00", "#FFF600", "#FFFF02", "#FFFF12", 
            "#FFFF22", "#FFFF32", "#FFFF42", "#FFFF52", "#FFFF62", 
            "#FFFF72", "#FFFF81", "#FFFF91", "#FFFFA1", "#FFFFB1", 
            "#FFFFC1", "#FFFFD1", "#FFFFE1", "#FFFFF1")
  if(n==1){ret<-orig[1];return(ret)}
  if(n==2){ret<-c(orig[1],orig[64]);return(ret)}
  rgb.hot <- t(col2rgb(orig))
  out<-apply(X = rgb.hot,2,function(a){round(spline(x = 1:64,y = a,n = n,method= "natural")$y)})
  out[out<0]<-0
  out[out>255]<-255
  ret<-rgb(out[, 1], out[, 2], out[, 3], maxColorValue = 255)
  return(ret)
}
