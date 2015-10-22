\name{viewR}
\alias{viewR}
\title{
Interactive Medical Image Viewer
}
\description{
This function allows for interactive viewing of 3D and 4D medical Images.
}

\usage{
viewR(data=NULL,overlay=NULL,otherData=NULL,xyz=NULL,ret=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{data}{
A 3D/4D array. If the array has an attribute named pixdim viewR will try and adjust the aspect appropriately. If NULL an interactive file selector will open to select the images one would like to view. If file paths to 3D/4D images are supplied they will be read and displayed. Multiple 3D files will be concatenated if the dimensions match. viewR can then be ran in movie mode. 
}

\item{overlay}{
A nifti image/file path to nifti image of the same dimensions which will be overlayed.
}

\item{otherData}{
A nifti image of the same dimensions which will also be displayed.
}
\item{xyz}{
numeric vector indicating initial voxel coordinates.
}
\item{ret}{
Logical indicating if read data should be returned.
}
}
\value{
If ret is set to TRUE the data will be returned.
}
\details{
Currently this viewR displays "world" coordinates as well as voxel coordiantes. The "world" coordinates are obtained by transforming the voxel coordinates using the qform matrix in the nifti header. Note that overlays and registration checks are not automatically displayed. They need to be checked in their respective menus.
}
\examples{
 \dontrun{
file<-system.file("extdata","mni.nii.gz",package="FIACH")
viewR(file)
}
}