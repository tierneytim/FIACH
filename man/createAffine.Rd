\name{createAffine}
\alias{createAffine}
\title{
Create Affine Matrix
}
\description{
This function creates an Affine matrix. The order of the matrix construction is Rotation-->Skew-->Scale--Translate}
\usage{
createAffine(translation = c(0, 0, 0), scales = c(1, 1, 1), skews = c(0, 0,
  0), angles = c(0, 0, 0), source = NULL, target = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{translation}{Translations along x, y and z in mm}

\item{scales}{Scale factors in x,y and z}

\item{skews}{Skews in the XY, XZ and YZ planes.}

\item{angles}{Pitch(rotation about x axis), roll(rotation about y axis) and yaw(rotation about z axis) angles, in radians.}

\item{source}{The source image for the transformation. This image determines the dimensions for the transformed image. If target is unspecified.}

\item{target}{The target image for the transformation. If \code{NULL} (the
default), it will be equal to \code{source}, or a rescaled version of it
if any of the \code{scales} are not 1. In the latter case the scales will
be reset back to 1 to produce the right effect.}
}


\value{
A 4x4 affine matrix representing the  transformation. The transform is from source to target apace. Note the output dimensions and qforms are stored as attributes of the matrix. These are used to construct the appied transformation matrix which is created in the following way: \code{solve(targetQform)\%*\%affine\%*\%sourceQform} 
}

\examples{
file<-system.file("extdata","mni.nii.gz",package="FIACH")
source <- readNii(file)

translation =c(6.4,-3.1,9.34)
angles = c(-.01,.3,.1)
skews = c(.1,.03,-.2)
scales = c(1.13,1.12,1.17)


affReg<-createAffine(source = source,
                    translation =translation,
                    angles = angles,
                    skews = skews,
                    scales = scales)

}
