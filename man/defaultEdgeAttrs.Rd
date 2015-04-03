\name{makeDefautEdgeData}
\alias{makeDefaultEdgeData}

\title{
Creates auxiliary data needed for SPIA method
}
\description{
This function creates a list containing auxiliary data needed in SPIA method for conversion between edge types and dividing interaction into three categories: positive, negative and neutral
}
\usage{makeDefaultEdgeData}
\value{  
A list of two data frames explained in the \emph{Details}
   \code{
  The format is:
List of 2
 $ graphite2SPIA: chr [1:26, 1:2] "binding" "control(In(ACTIVATION))" "control(In(INHIBITION))" "control(Out(ACTIVATION))" ...
  ..- attr(*, "dimnames")=List of 2
  .. ..$ : NULL
  .. ..$ : chr [1:2] "type" "spiaType"
 $ beta         :'data.frame':	25 obs. of  2 variables:
  ..$ rel : chr [1:25] "activation" "compound" "binding/association" "expression" ...
  ..$ beta: num [1:25] 1 0 0 1 -1 1 0 -1 -1 0 ...
} }
\details{
The first slot called \code{graphite2SPIA} contains a mapping table between edge types in topologies from \code{graphite} and edge types which are used in the implementation of SPIA in \code{SPIA} package. All of the edge types present in the topologies must be also covered by this table otherwise the method could not be applied. 


The second slot called \code{beta} divides the 25 interaction types into three categories: positive (\code{beta=1}), negative (\code{beta=-1} and neutral (\code{beta=0}) in the sense of gene regulation. Only user familiar with all the details of SPIA should change this. 
}
\source{
The data are manualy cerated from the unexported objects from \code{graphite} package version 1.10.1.
}


\examples{
str(makeDefaultEdgeData())
}
\keyword{manip}
