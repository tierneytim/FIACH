\name{viewNew}
\alias{viewNew}
\title{
Interactive Medical Image Viewer
}
\description{
This function allows for interactive viewing of 3D and 4D medical Images.
}

\usage{
viewNew(data=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{data}{
A 3D/4D array. If the array has an attribute named pixdim viewR will try and adjust the aspect appropriately. If NULL an interactive file selector will open to select the images one would like to view. If file paths to 3D/4D images are supplied they will be read and displayed. Multiple 3D files will be concatenated if the dimensions match. viewR can then be ran in movie mode. 
}
}

\value{
Nothing is returned.
}
\details{
Currently this viewR displays "world" coordinates as well as voxel coordiantes. The "world" coordinates are obtained by transforming the voxel coordinates using the qform matrix in the nifti header. 
}
\examples{
 \dontrun{
file<-system.file("extdata","mni.nii.gz",package="FIACH")
viewNew(file)
}
}