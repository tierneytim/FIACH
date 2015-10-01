\name{t2Grey}
\alias{t2Grey}

\title{
Transverse Relaxation
}
\description{
This function calculates the T2 or R2 of grey matter at varying field strength
}
\usage{
t2Grey(B0,relax)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{B0}{
Magnetic field strength specified in Tesla.
}

  \item{relax}{
Logical indicating if R2 or T2 should be returned
}

}

\value{
Either the T2 in ms or the R2 in Hz is returned
}
\author{
Tim Tierney
}
\examples{
t2Grey(1.5,relax=TRUE)
t2Grey(1.5,relax=FALSE)
}
