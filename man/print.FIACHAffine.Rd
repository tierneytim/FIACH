\name{print.FIACHAffine}
\alias{print.FIACHAffine}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Prints FIACHAffine object.}

\description{
prints the affine matrix
}
\usage{
\method{print}{FIACHAffine}(x,...)
}

\arguments{
  \item{x}{
An object of class FIACHAffine.
}
\item{...}{
Ignored
}
}
\author{
Tim Tierney
}

\examples{
file<-system.file("extdata","mni.nii.gz",package="FIACH")
source <- readNii(file)
aff<-createAffine(source=source)
print(aff)
}
