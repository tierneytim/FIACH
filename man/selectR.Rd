\name{selectR}
\alias{selectR}
\title{
Interactive File Selector
}
\description{
This function allows for interactive file selection using regular expressions.
}

\usage{
selectR(update=TRUE)
}
\arguments{
\item{update}{logical indicating whether or not the working directory should be changed to the current selected directory. This means that if this function is called a second time it will initiate in the directory it was in when previously called.
}
}

\value{
Returns a character vector of file paths. If no files are selected this will still be a character vector but will have length = 0. This function does not work well with network drives. If network drives are found this function initiates at the drive select section.
}

\examples{
 \dontrun{
selectR()
}
}