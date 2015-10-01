\name{getDatatype}
\alias{getDatatype}

\title{
Determine Datatype
}
\description{
This function returns the largest datatype from a set of files or images. The return value is a code in RNiftyReg format or the exact nifti code.
}
\usage{
getDatatype(input,type="RNiftyReg")
}
\arguments{
  \item{input}{
A character vector of file paths, an image or a list of images.
}
  \item{type}{
  A character vector. If equal to "RNiftyReg" the code will be returned in RNiftyReg format. If anything else the nifti code will be returned.
}

}
\value{
RNiftyReg datatype code or nifti datatype code.
}
\author{
Tim Tierney
}

\examples{
file<-system.file("extdata","motion_ex.nii.gz",package="FIACH")
getDatatype(file)
}