\name{hotMetal}
\alias{hotMetal}
\title{
Hot Metal
}
\description{
This function creates a hot metal colour palette}
\usage{
hotMetal(X)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
n must be an integer.
}
}
\value{
returns the colour palette
}

\examples{
plot.new()
pal<-as.matrix(hotMetal(200))
rasterImage(pal,0,0,1,1)  
}
