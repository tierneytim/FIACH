\name{gmm}
\alias{gmm}
\title{Gaussian Mixture Model}
\description{
EM algorithm for mixtures of gaussians}
\usage{
gmm(x, k, imeans = NULL, isd = NULL, ilambda = NULL, print=FALSE, tol=1e-8, maxit = 1000L)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{x}{numeric vector containing data}

\item{k}{integer specifying number of gaussians to be fitted}

\item{imeans}{optional initial means for mixture model}

\item{isd}{optional initial standard deviations for mixture model}

\item{ilambda}{optional initial mixing proportions for mixture model. Must sum to 1.}

\item{print}{logical indicating whether or not progress is printed to screen.}

\item{tol}{Convergence criteria for mixture model.}

\item{maxit}{integer specifying maximum number of iterations allowed}
}


\value{
returns a mixture model.
}

\examples{
test<-c(rnorm(1000),rnorm(1000,mean = 3,sd = 1))
hist(test)
a<-gmm(test,2)
plot(a)

}