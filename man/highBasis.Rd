\name{highBasis}
\alias{highBasis}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
High pass filter basis set
}
\description{
Generates cosines that can be used to high-pass filter data in a regression framework.
}
\usage{
highBasis(N, freq, tr)
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
ts.plot(highBasis(300,128,2.16))

}