\name{reslice}
\alias{reslice}
\title{
Reslice Image}
\description{
This function creates a 2D matrix from a  4D/3D image array.}
\usage{
reslice(data,newDim,currentDim=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
data must be a 4d or 3d array}
  \item{newDim}{
A length 3 numeric vector indicating the pixel dimensions of the resliced array}
  \item{currentDim}{
The current pixel dimensions. Only used if the array is not a nifti image}
}
\value{
returns a resliced array
}

\examples{
data<-readNii(system.file("extdata","mni.nii.gz",package="FIACH"))
re<-reslice(data,newDim = c(3,3,3))
}
