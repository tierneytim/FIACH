\name{rowsd}
\alias{rowsd}
\title{
Compute row standard deviations
}
\description{
This function computes the standard deviation of each row in a matrix.}
\usage{
rowsd(X)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
X must be a numeric matrix.
}
}
\value{
Returns a  vector containing the row standard deviations.
}

\examples{
\dontrun{
mat<-matrix(rnorm(100*100),ncol=100)
a<-rowsd(mat)


}
}