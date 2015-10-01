\name{colMad}
\alias{colMad}
\title{
Compute Column Median Absolute Deviations
}
\description{
This function computes the median absolute deviation of each column in a matrix.}
\usage{
colMad(X)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
X must be a numeric matrix.
}
}
\value{
returns a vector containing the column Median Absolute Deviations.
}

\examples{

mat<-matrix(rnorm(100*100),ncol=100)
a<-colMad(mat)

}
