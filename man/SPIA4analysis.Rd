\name{SPIA4analysis}
\alias{SPIA4analysis}
\alias{SPIA2graph}

\title{
Function to further process pathways for SPIA
}
\description{
The functions furthee process pathways prepared by \code{parseKEGG()} by \code{"SPIA"} method. \code{SPIA4analysis} prepares pathways for \code{SPIA()} function. \code{SPIA2graph} creates graphs, that can be used to visualize the pathway structure. Note, that the infomation about the edge type is lost during creation. 
}
\usage{
SPIA4analysis(pathways, EdgeAttrs = NULL)
SPIA2graph(pathways)
}

\arguments{
  \item{pathways}{
A list of pathways as obtained from \code{parseKEGG}
}
  \item{EdgeAttrs}{
A list of edge attributes. If \code{NULL}, it is created by \code{makeDefaultEdgeData()}. 
}
}
\value{
List of pathways
}
\author{
Ivana Ihnatova
}

\seealso{
\code{\link{parseKEGG}}
}
\examples{
\dontrun{
spiapt<-parseKEGG("./", "SPIA")
spiapt<-SPIA4analysis(spiapt)
SPIA(GEDM, group, spiapt, type="MA", preparePaths=FALSE)
}}
\keyword{ manip }

