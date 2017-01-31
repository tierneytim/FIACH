\name{applyAffine}
\alias{applyAffine}
\title{
Apply an Affine Matrix
}
\description{
This function Transforms the \code{source} image using a prespecified affine matrix.}
\usage{
applyAffine(source,aff)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{source}{The source image for the transformation.}
\item{aff}{an Affine matrix as porduced by \code{\link{createAffine}}  }
}


\value{
A numeric array containing the transformed image is returned
}

\examples{
source <- readNii(file)

translation =c(6.4,-3.1,9.34)
angles = c(-.01,.3,.1)
skews = c(.1,.03,-.2)
scales = c(1.13,1.12,1.17)


affReg<-createAffine(source = source,
                    translation =translation,
                    angles = angles,
                    skews = skews,
                    scales = scales)
                    
reg<-applyAffine(source,affReg)
invAffReg<-inverseAffine(affReg)
ireg<-applyAffine(reg,invAffReg)

diff<-ireg-source
sd(diff)
image(diff[,,45],col=grey(1:100/100),useRaster=TRUE)


}
