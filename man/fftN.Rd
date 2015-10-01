\name{fftN}
\alias{fftN}
\title{
Zero Padded 1D Fourier transform
}
\description{
This function is a simple wrapper of Armadillo's fft function. It allows for fast and easy zero padding of a signal. 
}
\usage{
fftN(X,N=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
X a numeric vector or matrix
}
  \item{N}{
Length of zero padded signal. If NULL the function will automatically pad sufficiently for a fast transform.
}
}
\value{
returns the Fourier transform of the signal.
}

\examples{
x<-matrix(rnorm(101*1000),nrow = 101,ncol = 1000)
system.time(test1<-fftN(x))

}