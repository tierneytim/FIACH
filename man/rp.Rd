\name{rp}
\alias{rp}
\docType{data}
\title{
Realignment parameters
}
\description{
Realignment Parameters which do no have any obvious extreme movement despite severe image artefact.
}
\usage{data(rp)}
\format{
  A data frame with 104 observations on the following 6 variables.
  \describe{
    \item{\code{V1}}{X Translation parameter}
    \item{\code{V2}}{Y Translation parameter}
    \item{\code{V3}}{Z Translation parameter}
    \item{\code{V4}}{Pitch}
    \item{\code{V5}}{Roll}
    \item{\code{V6}}{Yaw}
  }
}

\examples{
data(rp)
z<-rp[,3]                                                     
dz<- c(0,diff(z))                                                
plot(1:104,dz,ylim=c(-1,1),main="Volume-Volume Displacement",   
xlab="Time(scans)",ylab="Displacement(mm)")
lines(1:104,dz,col="blue",lwd=3)
legend("topright",legend="dZ/dt Translation",lwd=3,col="blue")
arrows(54,.5,91,-.184)
text(48,.5,"-.184mm")
}
\keyword{datasets}
