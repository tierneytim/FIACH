\name{spmOrth}
\alias{spmOrth}
\title{
Serial Orthogonalisation 
}
\description{
This function recursively orthogonalises variables in a matrix to the first column. It is designed to emulate the SPM approach to orthogonalisation. However instead of computing the results using a pseudo-inverse R's qr functions are used.
}
\usage{
spmOrth(a)
}
\arguments{
  \item{a}{
matrix to be orthogonalised
}
}

\value{
The recursively orthogonalised matrix is returned
}
\author{
Tim Tierney
}


\examples{
a<-basisFunctions(1/16,orth=FALSE)
ts.plot(a,ylim=c(-.007,.015))
ts.plot(spmOrth(a),ylim=c(-.007,.015))
}
