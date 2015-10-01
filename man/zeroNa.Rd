\name{zeroNa}
\alias{zeroNa}
\title{
Zero non finite values
}
\description{
This function replaces NA, Inf, -Inf, and NaN values of a matrix, array or vector with zero.
}

\usage{
zeroNa(input)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{input}{
A numeric vector, matrix, or array. If an array it must have less than 5 dimensions. 
}
}
\value{
Returns the same input but with zeros in place of the non-finite values.
}

\examples{
 test<-matrix(rnorm(100),10,10)

test[5,5]<-NA
test
zeroNa(test)

test[5,5]<-NaN
test
zeroNa(test)

test[5,5]<- -Inf
test
zeroNa(test)

test[5,5]<- +Inf
test
zeroNa(test)
}