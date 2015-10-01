\name{sinc}
\alias{sinc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Sinc Pulse}
\description{
Creates a sinc with desired properties.
}
\usage{
sinc(fh,fl,tw,sf,type,n)
}
\arguments{
  \item{fh}{
High pass frequency in Hz.
}
  \item{fl}{
  Low pass frequency in Hz.
}
\item{tw}{
Transition width in Hz.
}

\item{sf}{
sampling frequency in Hz.}

\item{type}{
string of either "low","high","band" depending on what you wish to do.}
\item{n}{
range over which filter is to be evaluated.}
}

\value{
outputs a sinc pulse evaluated over n.}

\author{
Tim Tierney
}

\examples{
par(mfrow=c(2,1))
si<-sinc(fl=.08,tw=0,type="low",n=-40:40,sf=1/2.16)
ts.plot(si,ylab="",main="Sinc Pulse")
s<-rawPeriodogram(si,1/2.16)
}