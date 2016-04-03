kmeansMask<-function(x,logData=FALSE,retFit=FALSE){
  if(logData){
    x[x<=0]<-1
    x<-log(x)
    }
  resamp<-quantile(x,probs = seq(0,1,length.out = 1000))
  mod<-kmeans(resamp,2)
  thresh<-resamp[which.max(abs(diff(mod$cluster)))]
  mask<-ifelse(x<=thresh,0,1)
  if(retFit){mask<-list(mask=mask,fit=mod$betweenss/mod$totss)}
  return(mask)
}


