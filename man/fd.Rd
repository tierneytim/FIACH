\name{fd}
\alias{fd}
\title{
Framewise Displacement
}
\description{
This function computes the Framewise Displacement of Realignment Parameters. This function assumes that the first 6 columns of the input are the volumetric realignment parameters.}
\usage{
fd(input)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{input}{
A character vector specifying the path to the rp file or a matrix.
}
}
\value{
returns a time series containing the fd across time
}

\examples{
data(rp)
fd(rp)
}
