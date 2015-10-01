\name{sepConvolve3d}
\alias{sepConvolve3d}
\title{
Fast 3d Convolution
}
\description{
This function  convolves data in three different directions with potentially different impulses. 
}

\usage{
sepConvolve3d(x,kernX, kernY, kernZ)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A 3D/4D array. If a 4D array then each 3d array along the fourth dimension is convolved.
}
  \item{kernX}{
filter kernel to be used in the x direction
}
  \item{kernY}{
filter kernel to be used in the y direction
}
 \item{kernZ}{
filter kernel to be used in the z direction
}
}
\value{
returns the convolved array
}

\examples{
func<-readNii(system.file("extdata","motion_ex.nii.gz",package="FIACH"))

kernx<-gaussKernel(8,21,3.33)
kerny<-gaussKernel(8,21,3.33)
kernz<-gaussKernel(8,21,4)

smooth<-sepConvolve3d(func,kernx,kerny,kernz)
}
