\name{lowBasis}
\alias{lowBasis}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Low pass filter basis set
}
\description{
Generates cosines that can be used to low-pass filter data in a regression framework.
}
\usage{
lowBasis(N, freq, tr)
}

\arguments{
  \item{N}{
Number of time-points to be filtered.
}
  \item{freq}{
This is the cutoff frequency. It is specified in seconds for comparability with SPM.
}
  \item{tr}{
This is the time between samples(or the TR for fMRI data.)
}
}

\value{
A matrix containing the basis set.
}

\author{
Tim Tierney
}

\examples{

## 300 scans, 10 second cutoff, tr=2.16 seconds
lb<-lowBasis(300,10,2.16)
dim(lb)
plot.ts(lb[,1:10])

}