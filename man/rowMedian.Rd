\name{rowMedian}
\alias{rowMedian}
\title{
Compute Row Medians
}
\description{
This function computes the median of each row in a matrix.}
\usage{
rowMedian(X)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
A must be a numeric matrix.
}
}
\value{
Returns a vector containing the medians.
}

\examples{

mat<-matrix(rnorm(100*100),ncol=100)
a<-rowMedian(mat)
}