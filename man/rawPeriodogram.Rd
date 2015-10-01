\name{rawPeriodogram}
\alias{rawPeriodogram}
\title{
Raw Periodogram
}
\description{
This function computes and plots a raw periodogram. This function may produce inaccurate results on stochastic data. 
}
\usage{
rawPeriodogram(x,sf=NULL,plot=TRUE,amp=FALSE,phase=FALSE,N=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
x a numeric vector
}
\item{sf}{
Sampling frequency in Hz. If specified the power spectrum will be plotted with Hz instead of normalised frequency on the x-axis. 
}
\item{plot}{
Whether or not to plot the periodogram. If more than one series is used only the first will be plotted.
}
\item{amp}{
Whether or not to return the amplitude spectrum.
}

\item{phase}{
Whether or not to return the phase spectrum.
}
\item{N}{
Length of Fourier transform used to compute spectrum. If NULL N will be chosen for speed and not to prevent scalloping loss.
}

}

\value{
Returns the periodogram and optionally plots it.
}

\examples{

hz.slow<-5
hz.fast<-50
t<-seq(0,1,length.out=300)

sin.sig.slow<-sin(2*pi*t*hz.slow)
sin.sig.fast<-sin(2*pi*t*hz.fast)
sin.sig.combo<-sin.sig.slow+sin.sig.fast

ts.plot(sin.sig.combo)
rawPeriodogram(sin.sig.combo)
rawPeriodogram(sin.sig.combo,300)

}
