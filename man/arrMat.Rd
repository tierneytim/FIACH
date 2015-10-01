\name{arrMat}
\alias{arrMat}
\title{
Array to matrix transformation.
}
\description{
This function creates a 2D matrix from a  4D/3D image array.}
\usage{
arrMat(input)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{input}{
input must be a 4d or 3d array}
}
\value{
returns a matrix(time x No. of voxels)
}

\examples{
arr<-array(dim=c(64,64,30,10)) ## 4D array
mat<-arrMat(arr)
dim(mat)

}
