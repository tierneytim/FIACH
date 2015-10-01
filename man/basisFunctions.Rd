\name{basisFunctions}
\alias{basisFunctions}
\title{
Informed basis set
}
\description{
This function produces the impulse response for the canonical haemodynamic response function and its derivatives in time and space at a desired sampling rate.
}
\usage{
basisFunctions(RT,orth=TRUE)
}
\arguments{
  \item{RT}{
RT is the time between samples. It must be specified in seconds (e.g. 1/16).
}
\item{orth}{
Boolean describing whether or not the basis should be orthogonalised.
}
}

\value{
A matrix is returned containing all three impulse responses.
}
\author{
Tim Tierney
}
\examples{
RT=1/16
basis<-basisFunctions(RT=RT)
basis<-basis/max(basis)
x<-seq(0,32,RT)
ylim<-c(min(basis)-sd(basis),max(basis)+sd(basis))
plot(x,basis[,1],lwd=7,col="red",type="l",
      main="Informed Basis Set",xlab="Time(seconds)",
      ylab="Intensity(A.U.)",ylim=ylim)
      
y<-seq(min(ylim),max(ylim),length.out=10)
abline(h=y,col="grey")
legend(x="topright",                                 
       legend=c("Canonical HRF",
                "Temporal Derivative",
                "Dispersion Derivative"),
       col=c("red","blue","green"),lwd=2)
lines(x,basis[,1],lwd=7,col="red")
lines(x,basis[,2],lwd=7,col="blue")
lines(x,basis[,3],lwd=7,col="green")

}
