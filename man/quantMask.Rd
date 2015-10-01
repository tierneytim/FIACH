\name{quantMask}
\alias{quantMask}
\title{
Compute the Mask Vector}
\description{
This function creates a mask vector for an image using a user specified quantile.}
\usage{
quantMask(x,quant=.7)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
x can be a 4D/3D array, 2D matrix or 1D vector representing an image.
}
  \item{quant}{
 The quantile level used for masking. Defaults to .7.
}
}
\value{
returns a binary vector of zeros and ones.
}

\examples{
arr<-array(rnorm(10*10*10*10),dim=c(10,10,10,10))
mask<-quantMask(arr,.7)

}