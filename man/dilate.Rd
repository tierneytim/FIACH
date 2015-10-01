\name{dilate}
\alias{dilate}
\title{
Separable Flat dilation
}
\description{
This function dilates data with a flat kernel in up to three dimensions. 
}

\usage{
dilate(input,k)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{input}{
A vector/matrix or 3D array. 
}
  
 \item{k}{
Window width. Note that the window is symmetric around the i^th element with a total of 2*k+1 elements.
}
}
\value{
returns the dilated array
}

\examples{
arr<-array(rnorm(30*30*30), dim=c(30,30,30))
d<-dilate(arr, k=1)

}
