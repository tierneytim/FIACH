\name{hampel}
\alias{hampel}
\title{
Hampel Filter
}
\description{
This function Hampel Filters each column of a matrix.}
\usage{
hampel(x,k=3,t0=3)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Numeric Matrix where the columns are to be filtered. 
}
  \item{k}{
Window width. Note that the window is symmetric around the i^{th} element with a total of 2*k+1 elements. 
}
  \item{t0}{
 Threshold, expressed in  Median Absolute Deviations. Note that setting the threshold to 0 creates a Median filter. 
}

}
\value{
Returns a matrix containing the filtered data
}

\examples{
mat<-matrix(rnorm(100*100),ncol=100)
a<-hampel(mat)

}