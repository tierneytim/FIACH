\name{fiachMerge}
\alias{fiachMerge}
\title{
Combine 3D .nii files into 4D .nii files
}
\description{
This function Combine 3D .nii files into 4D .nii files
}
\usage{
fiachMerge(input,compress=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{input}{
Character string specifying filepaths to 3D files.}

  \item{compress}{
Logical indicating if files should be compressed.}
}
\value{
Nothing is returned to R but the 4D file is written to the directory from which the 3D files came from.
}

\examples{

}
