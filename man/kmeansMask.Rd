\name{kmeansMask}
\alias{kmeansMask}
\title{
Binarize vector, matrix or array}
\description{
This function creates a mask for an image using kmeans clustering.}
\usage{
kmeansMask(x,logData=FALSE,retFit=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
x can be a vector,matrix or array.
}
  \item{logData}{
Log transform data before applying the kmeans clustering. Must be TRUE or FALSE. 
}
  \item{retFit}{
Determine whether the fit should be returned. Must be TRUE or FALSE.
}
}
\value{
returns a binary vector, matrix or array of zeros and ones.
}

\examples{
file<-system.file("extdata","motion_ex.nii.gz",package="FIACH")
arr<-readNii(file)[,,,1]
mask<-kmeansMask(arr)

}