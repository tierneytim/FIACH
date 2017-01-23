\name{getAnalyzeOrigin}
\alias{getAnalyzeOrigin}
\title{
Analyze origin
}
\description{
This function returns the image origin having being supplied a file path}
\usage{
getAnalyzeOrigin(fname)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fname}{
Path to the analyze file. The path can either be to the .img file .hdr file. Compression using gzip is supported. 
}
}
\value{
Returns a  vector of length three containing the origin. This may be nonsensical if the origin has not been set. This may be the case in a number of files,
}

\examples{
file<-system.file("extdata","avg152T1.hdr.gz",package="FIACH")
getAnalyzeOrigin(file)

}
