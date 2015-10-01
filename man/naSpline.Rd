\name{naSpline}
\alias{naSpline}
\title{
Interpolate NAs
}
\description{
This function interpolates NAs in matrices along columns.}
\usage{
naSpline(mat,maxgap=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mat}{
A must be a numeric matrix.
}
  \item{maxgap}{
positive integer indicating over what gaps values can be interpolated. Longer gaps will be less accurate.
}
}
\value{
Returns the interpolated matrix.
}

\examples{
x<-rnorm(100)
x[20]<-30
naData<-badData(x)
fixed<-naSpline(naData)
par(mfrow=c(2,1))
ts.plot(x)
ts.plot(fixed)

}
