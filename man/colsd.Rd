\name{colsd}
\alias{colsd}
\title{
Compute column standard deviations
}
\description{
This function computes the standard deviation of each column in a matrix.}
\usage{
colsd(X)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
X must be a numeric matrix.
}
}
\value{
Returns a  vector containing the column standard deviations.
}

\examples{
mat<-matrix(rnorm(100*100),ncol=100)
a<-colsd(mat)
}
