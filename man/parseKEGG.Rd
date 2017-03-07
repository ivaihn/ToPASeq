\name{parseKEGG}
\alias{parseKEGG}

\title{
Function to covert KEGG pathways for topology-based pathway analysis
}
\description{
The function processes KGML files from KEGG located on local drive to the format suitable for topology-based pathway analysis. The processing is method-specific and the processed object can be used only with the matching method.}
\usage{
parseKEGG(pathfolder, method)
}
\arguments{
  \item{pathfolder}{
A character, path to the folder with KGML files. For working directory use \code{"./"}
}
  \item{method}{
A character, name of the method. Currently implemented are: \code{SPIA}, \code{PRS}, \code{CePa}, \code{DEGraph}, \code{clipper}
}
}
\details{
Nodes, relations and reactions found in KGML files are processed as in the original implementation of the methods, or as stated in the original publication. For clipper, the \code{graphite}'s approach described in Sales, 2012 was used. 
}
\warning{
\code{preparePaths=FALSE} must be used in functions that apply the method e.g. \code{SPIA()}
}
\value{
List of parsed pathways. The exact form is method-specific. 
}
\author{
Ivana Ihnatova
}
\examples{
\dontrun{
# to parse all xml files in the working directory
parseKEGG("./", "PRS")
}
}

\keyword{ manip }% use one of  RShowDoc("KEYWORDS")

