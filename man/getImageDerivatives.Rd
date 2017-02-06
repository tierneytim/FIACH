\name{getImageDerivatives}
\alias{getImageDerivatives}
\title{
Image derivatives
}
\description{
This function calcualtes image derivatives for 3d images}
\usage{
getImageDerivatives(im,translation=TRUE,rotation=TRUE,skews=FALSE,scales=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{im}{
im must be a 3d array or 3d niftiImage.
}
  \item{translation}{
Boolean indicating whether the derivative with respect to translation should be returned.
}
  \item{rotation}{
Boolean indicating whether the derivative with respect to rotation should be returned.
}
  \item{skews}{
Boolean indicating whether the derivative with respect to skew should be returned.
}
  \item{scales}{
Boolean indicating whether the derivative with respect to scale should be returned.
}
}

\value{
Returns a list containing the derivative images.
}

\examples{
file<-system.file("extdata","motion_ex.nii.gz",package="FIACH")
im<-readNii(file)
Deriv<-getImageDerivatives(im[,,,1])
}
