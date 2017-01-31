\name{inverseAffine}
\alias{inverseAffine}
\title{
Create  inverse of Affine Matrix
}
\description{
This function inverts an Affine matrix. }
\usage{
inverseAffine(aff)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{aff}{an Affine matrix as porduced by \code{\link{createAffine}}  }
}


\value{
A 4x4 affine matrix representing the  transformation. The transform is from source to target apace. Note the output dimensions and qforms are stored as attributes of the matrix. These are used to construct the appied transformation matrix which is created in the following way: \code{solve(targetQform)\%*\%affine\%*\%sourceQform} 
}

\examples{
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
                    
invAffReg<-inverseAffine(affReg)
affReg\%*\%invAffReg
}
