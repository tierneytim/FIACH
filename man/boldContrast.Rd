\name{boldContrast}
\alias{boldContrast}

\title{
Maximum Theoretical BOLD signal observable
}
\description{
This function returns a theoretical threshold beyond which BOLD signal changes are likely to be artefactual based on the models presented in Yablonskiy and Haacke(1994).
}
\usage{
boldContrast(B0, te, plot = TRUE, random = TRUE, alpha = 0.38, hct = 0.4, cbf.base = 55, 
             chi = 1.8e-07, e.base = 0.4, w0 = 267500000,e.act=.1,cbf.act=2*cbf.base)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{B0}{
Magnetic field strength specified in Tesla.
}
  \item{te}{
Echo Time  specified in ms
}
  \item{plot}{
Logical indicating whether or not to plot results.
}
  \item{random}{
Logical indicating if the theoretical model assumes the cylinders embedded in the medium are randomly orientated. If false the cylinders are assumed to be in parallel.
}
  \item{alpha}{
Power to which cbf.base is raised to estimate Cerebral Blood Volume. }
  \item{hct}{
Hematocrit expressed as a real number between 0 and 1. 
}
  \item{cbf.base}{
Cerebral Blood Flow at rest, specified in ml per mg per minute.
}
  \item{chi}{
Volume susceptibility difference between oxygenated and de-oxygenated blood specified in cgs units.
}
  \item{e.base}{
  1- oxygenation fraction during rest. ranging between 0 and 1
}
  \item{w0}{
Gyromagnetic ratio of the proton expressed in radians per second per tesla
}
 \item{e.act}{
  1- oxygenation fraction during activation. ranging between 0 and 1. Defaults to .1 to produce a conservative threshold. A value of .2 might produce a more realistic threshold.
}
 \item{cbf.act}{
  Cerebral Blood Flow during activation, specified in ml per mg per minute
}
}

\value{
A threshold expressed in percent signal change is returned 
}
\author{
Tim Tierney
}
\examples{
## magnetic field strength of 1.5T and TE of 30ms
## assuming parallel orientation of cylinders.
boldContrast(1.5,30,plot=TRUE,random=FALSE) 

## magnetic field strength of 1.5T and TE of 
## 30ms assuming random orientation of cylinders
boldContrast(1.5,30,plot=TRUE,random=TRUE)   

}
