\name{rowMad}
\alias{rowMad}
\title{
Compute row Median Absolute Deviations
}
\description{
This function computes the Median Absolute Deviation of each row in a matrix.}
\usage{
rowMad(X)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
X must be a numeric matrix.
}
}
\value{
Returns a  vector containing the row Median Absolute Deviations.
}

\examples{
mat<-matrix(rnorm(100*100),ncol=100)
a<-rowMad(mat)

}
