\name{plot.gmm}
\alias{plot.gmm}
\title{
Plot a Gaussian Mixture Model
}

\description{
Plots an object of class gmm(Gaussian Mixture Model).
}
\usage{
\method{plot}{gmm}(x,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
an object of class gmm
}
 \item{...}{Graphical parameters passed to \code{plot} command.}
}
\examples{

test<-c(rnorm(1000),rnorm(1000,mean = 3,sd = 1))
hist(test)
a<-gmm(test,2)
plot(a)

}