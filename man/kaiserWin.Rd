\name{kaiserWin}
\alias{kaiserWin}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Kaiser Window}
\description{
Creates a windowed sinc pulse using a kaiser window.
}
\usage{
kaiserWin(fh,fl,tw,sf,d.sa,d.pbr,type)
}

\arguments{
  \item{fh}{
High pass frequency in Hz.}
  \item{fl}{
  Low pass frequency in Hz}
\item{tw}{
Transition width in Hz.}
\item{sf}{
Sampling frequency in Hz.}
\item{d.sa}{
Desired Stopband Attenuation specified in dB.}
\item{d.pbr}{
Desired Passband ripple specified in dB.}
\item{type}{
string of either "low","high","band" depending on what you wish to do.}
}

\value{
Outputs FIR filter coefficients to desired specifications.}

\author{
Tim Tierney
}

\examples{
kais<-kaiserWin(fl=.08,tw=.025,sf=1/2.16,d.sa=70,d.pbr=.1,type="low")
par(mfrow=c(2,1))
ts.plot(kais,main="Low Pass Kaiser Window")
rawPeriodogram(kais,1/2.16)
}
