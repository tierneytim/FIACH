\name{highPass}
\alias{highPass}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
SPM high pass filter
}
\description{
Removes low frequency drifts from data.
}
\usage{
highPass(x, freq, tr)
}

\arguments{
  \item{x}{
Numeric matrix or vector. If a vector is supplied it will be coerced to a matrix.
}
  \item{freq}{
This is the cutoff frequency. It is specified in seconds for comparability with SPM.
}
  \item{tr}{
This is the time between samples(or the TR for fMRI data.)
}
}

\value{
A high pass filtered matrix will always be returned.
}

\author{
Tim Tierney
}

\examples{
t<-seq(1,300)
y<-cos(.01*pi*t)+.1*cos(.2*pi*t)   ## time series with low frequency drift
plot(t,y,type="l",col="red",lwd=5) ## plot it
hp.y<-highPass(y,128,2.16)        ## filter it at 128s (like in SPM) with a tr of 2.16
lines(t,hp.y,col="blue",lwd=5)     ## result has no drift
legend(x="topright",               ## Create the Legend
legend=c("Raw Signal",
         "Highpass Filtered"),
col=c("red","blue"),lwd=2)

}