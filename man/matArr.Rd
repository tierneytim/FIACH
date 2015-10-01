\name{matArr}
\alias{matArr}
\title{
Matrix to Array Transformation.
}
\description{
This function creates a 4D/3D image array from a 2D matrix  }
\usage{
matArr(mat,dim)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mat}{
input must be a 2d matrix}

  \item{dim}{
dimensions of array}
}
\value{
returns an array of dimension dim
}

\examples{

mat<-matrix(nrow=50,ncol=30*30*15)  
arr<-matArr(mat,dim=c(30,30,15,50))
dim(arr)

}
