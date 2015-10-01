\name{readNii}
\alias{readNii}

\title{
Input function for fmri data
}
\description{
This function reads 4D and 3d nifti and nifti files
}
\usage{
readNii(input,fourD=TRUE)
}
\arguments{
  \item{input}{
A character vector containing a single string or multiple strings. If multiple strings are provided the data is concatenated and is assumed to be from the same session. The data can be in 4d.nii format or .img and .hdr pairs. If input input is in .img/.hdr pair only the .img file need be specified
}
  \item{fourD}{
If set to false nii objects are returned for each image.
}

}
\value{

A 4D array is returned with time being the fourth dimension
}
\author{
Tim Tierney
}

\examples{
file<-system.file("extdata","motion_ex.nii.gz",package="FIACH")
func<-readNii(file)

}