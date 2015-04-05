\name{KEGG2pathway}
\alias{KEGG2pathway}

\title{
Function to parse KEGG KGML file into a pathway
}
\description{
The function parses a KGML file from KEGG into a pathway. 
}
\usage{
KEGG2pathway(file, expandGenes = TRUE, expandCom = TRUE, nongene = c("keep", "propagate", "discard"), ident = "KEGGnative", database = "KEGG", species = NULL)
}

\arguments{
  \item{file}{
Character, the name of the file to be parsed. Download manually or in bulk from KEGG
}
  \item{expandGenes}{
Logical, should multi-gene nodes be expanded into separate nodes?
}
  \item{expandCom}{
Logical, should undirected binding interactions be added between nodes from one group (usually multi-subunit protein complex, which is turned into a clique)
}
  \item{nongene}{
Character, how should be the non-gene nodes parsed? If \code{"discard"} they are removed from the pathway. If \code{"propagate"}, they are removed but the interactions are preserved (e.g. if gene A interacts with compound c and compound c interacts with gene B, then the interaction between A and B is preserved. Otherwise, they are kept in the pathway topology
}
  \item{ident}{
Character, the type of the node identifiers. 
}
  \item{database}{
Character, the name of the database
}
  \item{species}{
Character, the three-letter code for the species-specific pathways. If \code{NULL} then, the first 3 letters from the \code{file} are used.
}
}
\value{
A \code{pathway} }
\author{
Ivana Ihnatova
}


\seealso{
 \code{\link[graphite]{pathway-class}}
}
\keyword{ manip }

