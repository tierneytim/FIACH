\name{realign}
\alias{realign}

\title{
Realignment of functional MRI data.
}
\description{
This function realigns a series of 3d functional images to the first in the series.
}
\usage{
realign(files,quality=.9,write=TRUE,plotConvergence=FALSE)
}
\arguments{
  \item{files}{
A character vector specifying the paths to the files in need of realignment
}
  \item{quality}{
  A real number between 0 and 1 that determines how many voxels are used inthe registration. A value of 1 uses all the voxels.
}
  \item{write}{
  Logical indicating if the registered files should be written to the file system
}
  \item{plotConvergence}{
  Logical indicating if the convergence of the method used to select the most interesting voxels should be plotted.
}


}
\value{
A matrix containing the  estimated translations and rotations. 
}
\author{
Tim Tierney
}

\examples{
## TO DO
}