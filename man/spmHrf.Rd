\name{spmHrf}
\alias{spmHrf}
\title{
Canonical Haemodynamic Response Function.
}
\description{
This function produces the impulse response for the canonical haemodynamic response function at a desired sampling rate.
}
\usage{
spmHrf(RT, p = c(6, 16, 1, 1, 6, 0, 32))
}
\arguments{
  \item{RT}{
 RT is the desired sampling frequency. It must be specified in Hz(e.g. 1/16). 
}
  \item{p}{
  Default parameters for impulse response estimation identical to SPM's implementation.
}
}
\value{

  \item{comp1 }{The first component of the list is the impulse response.}
  \item{comp2 }{The second component are the parameters used to estimate it.}

}
\author{
Tim Tierney
}

\examples{
RT<-1/16
can.hrf<-spmHrf(RT)[[1]]
x<-seq(0,32,RT)
plot(x,can.hrf,lwd=7,col="red",type="l",
main="Canonical HRF",xlab="Time(seconds)",
ylab="Intensity(A.U.)")
}
