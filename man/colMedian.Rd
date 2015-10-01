\name{colMedian}
\alias{colMedian}
\title{
Compute Column Medians
}
\description{
This function computes the median of each column in a matrix.}
\usage{
colMedian(X)
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
a<-colMedian(mat)

}
