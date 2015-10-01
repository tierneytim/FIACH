\name{badData}
\alias{badData}
\title{
Identifies spurious observations 
}
\description{
This function identifies spurious observations using a mix of stochastic and deterministic thresholds.}
\usage{
badData(X, meds=NULL, mads=NULL, nMads=1.96, t=0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
X must be a numeric matrix.
}
  \item{meds}{
Vector containing column center. If NULL will calculate the median.
}
\item{mads}{
Vector containing column variability. If NULL will calculate the Median Absolute Deviation.
}
\item{nMads}{
The stochastic threshold specified in terms of MADs
}
\item{t}{
The deterministic threshold. Specified in percentage distance from median.
}
}
\value{
Returns a matrix with NaN where data was identified as spurious
}

\examples{
x<-rnorm(100)
x[20]<-30
badData(x)

}
