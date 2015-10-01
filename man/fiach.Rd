\name{fiach}
\alias{fiach}
\title{
FIACH 
}
\description{
This function makes use of The EM algorithm to segment and model noisy areas of the brain. It also utilizes a thresholded spline filter to clean spurious observations in the image.
}
\usage{
fiach(input, t, tr, rp=NULL, maxgap = 1, freq = 128,nMads = 1.96)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{input}{
A character vector containing a single string or multiple strings. If multiple strings are provided the data is concatenated and is assumed to be from the same session.The data can be in 4D.nii format or .img and .hdr pairs. If input input is in .img/.hdr pair only the .img file need be specified. These images should be realigned.
}
  \item{t}{
A threshold expressed in percent signal change beyond which the spline filter is used. This threshold is ideally computed using boldContrast.
}
  \item{tr}{
The time between scans. Needed for the high pass filtering. 
}
  \item{rp}{
The realignment parameter files that will be appended with additional regressors.  
}
  \item{maxgap}{
  max number of consecutive time-points you will allow to be filtered before scrubbing (replacement with median) takes place.
}
  \item{freq}{
Desired frequency to high pass filter data at. Specified in seconds. CAUTION!!! THE DEFAULT IS 128 WHICH MAY NOT BE APPROPRIATE FOR ALL DESIGNS.
}

  \item{nMads}{
Number of MADs used for noise calculation.
}

}


\value{
Nothing is returned to R but the filtered files are written to the directory they came from with the prefix filt_. A directory is also created to store the various diagnostic images and plots produced by this method(median, rTSNR, mask and rTSNR histogram). The regressors to be used in further analysis are in the noise_basis6 file. The global median signal is also outputted in gs.txt file. The framewise displacement is also outputted appended to the noise regressors in fd_noise.txt if movement regressors are supplied. 
}
\author{
Tim Tierney
}

\examples{
\dontrun{
#### create the necessary files ###

dir.create("fiach_example/")
file<-system.file("extdata","motion_ex.nii.gz",package="FIACH")
arr<-readNii(file)
RNiftyReg::writeNifti(arr,"fiach_example/motion_ex.nii.gz",datatype = "short")
data(rp)
write.table(rp[1:13,],"fiach_example/rp.txt",col.names=FALSE,row.names=FALSE)

# running FIACH #
t<-boldContrast(1.5,30)
tr<-2.16
system.time(fiach("fiach_example/motion_ex.nii.gz", t=t,tr=tr, rp="fiach_example/rp.txt"))

## Note that this is a toy example.
## The results are meaningless.
## This is only an example of how
## one would use fiach in a script.
}
}