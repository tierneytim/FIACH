\name{viewR}
\alias{viewR}
\title{
Interactive Medical Image Viewer
}
\description{
This function allows for interactive viewing of 3D and 4D medical Images.
}

\usage{
viewR(data=NULL,ret=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{data}{
A 3D/4D array. If the array has an attribute named pixdim viewR will try and adjust the aspect appropriately. If NULL an interactive file selector will open to select the images one would like to view. If file paths to 3D/4D images are supplied they will be read and displayed. Multiple 3D files will be concatenated if the dimensions match. viewR can then be ran in movie mode. Currently this viewR does not display world coordinates. Instead voxel coordinates are used. As such it does not enforce information obtained from nifti headers.
}
\item{ret}{
Logical indicating if read data should be returned.
}
}
\value{
If ret is set to TRUE the data will be returned.
}

\examples{
 \dontrun{
url <- "http://nifti.nimh.nih.gov/nifti-1/data/filtered_func_data.nii.gz"
urlfile <- tempfile(pattern="filtered_func_data", fileext=".nii.gz")
download.file(url, urlfile, quiet=TRUE)
data<-readNii(urlfile)
viewR(data)
## Try coordinate  47,25,10
}
}