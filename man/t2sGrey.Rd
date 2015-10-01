\name{t2sGrey}
\alias{t2sGrey}

\title{
Effective Transverse Relaxation}
\description{
This function returns a theoretical value for T2* based on the models presented in Yablonskiy and Haacke (1994).
}
\usage{
t2sGrey(B0, alpha = 0.38, hct = 0.4, cbf.base = 55, chi = 1.8e-07, 
        e.base = 0.4, w0 = 267500000,r2=t2Grey(B0,TRUE),relax=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{B0}{
Magnetic field strength specified in Tesla.
}

  \item{alpha}{
Power to which cbf.base is raised to estimate Cerebral Blood Volume. }
  \item{hct}{
Hematocrit expressed as a real number between 0 and 1. 
}
  \item{cbf.base}{
Cerebral Blood Flow at rest, specified in ml per mm per minute.
}
  \item{chi}{
Volume susceptibility difference between oxygenated and deoxygenated blood specified in cgs units.
}
  \item{e.base}{
  1- oxygenation fraction during rest. ranging between 0 and 1
}
  \item{w0}{
Gyromagnetic ratio of the proton expressed in radians per second per tesla
}
\item{r2}{Transverse relaxation rate. Specified in Hz}
\item{relax}{Boolean indicating whether R2* or T2* should be returned}
}

\value{
Either R2* in Hz or T2* in ms is returned}
\author{
Tim Tierney
}
\examples{
t2sGrey(3,relax=TRUE)
t2sGrey(3,relax=FALSE)

}
