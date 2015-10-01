\name{gaussKernel}
\alias{gaussKernel}
\title{
1 dimensional Gaussian smoothing kernel
}
\description{
This function produces a Gaussian smoothing kernel that can be used in separable convolution to smooth images.
}

\usage{
gaussKernel(fwhm=8,ksize=5,voxsize=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{fwhm}{Full width at half maximum of the kernel.}
\item{ksize}{Size of kernel. Must be Odd}
\item{voxsize}{Size of the voxel in the direction you wish to smooth in. If smoothing a non medical image leave voxsize as 1 unless anisotropic smoothing is required.}
}

\value{
Returns a Kernel
}

\examples{
ts.plot(gaussKernel(ksize=23,voxsize=1))
}
