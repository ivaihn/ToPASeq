\name{ToPASeq-package}
\alias{ToPASeq-package}
\alias{ToPASeq}
\docType{package}
\title{
Package for topology-based pathway analysis of microarray and RNASeq data
}
\description{
The package implementats  several methods for topology-based pathway analysis of microarray data. The methods present in here are: SPIA, TopologyGSA, DEGraph, Clipper, PWEA, TAPPA, TBS. SPIA, PWEA and TBS were also adapted for RNASeq data.
}
\details{
\tabular{ll}{
Package: \tab ToPASeq\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2014-03-04\cr
License: \tab AGPL-3\cr
}

}
\author{
Ivana Ihnatova

Maintainer: Ivana Ihnatova <ihnatova@iba.muni.cz>

}


\keyword{ package }

\examples{

\dontrun{
if (require(DEGraph)) {
  data("Loi2008_DEGraphVignette")
  pathways<-biocarta[1:10]
  SPIA(exprLoi2008, classLoi2008,pathways , type="MA", logFC.th=-1, IDs="entrez")
  DEGraph(exprLoi2008, classLoi2008, pathways, type="MA")
  TAPPA(exprLoi2008, classLoi2008, pathways, type="MA")
  TopologyGSA(exprLoi2008, classLoi2008, pathways, type="MA", nperm=200)
  Clipper( exprLoi2008, classLoi2008+1, pathways,type="MA", test="mean")
  PWEA(exprLoi2008, classLoi2008, pathways,  type="MA", nperm=100)
  TBS( exprLoi2008, classLoi2008, pathways, type="MA",  logFC.th=-1, nperm=100)
}
if (require(gageData)) {

 data(hnrnp.cnts)
 group<-c(rep("sample",4), rep("control",4))
 SPIA( hnrnp.cnts, group, biocarta[1:10], type="RNASeq",  logFC.th=-1, IDs="entrez", test="limma")
 DEGraph(hnrnp.cnts, group, biocarta[1:10], type="RNASeq", norm.method="TMM")
 TAPPA( hnrnp.cnts, group, biocarta[1:10], type="RNASeq", norm.method="TMM")
 TopologyGSA(hnrnp.cnts, group, biocarta[1:10], type="RNASeq",nperm=200, norm.method="TMM")
 Clipper(hnrnp.cnts, group,biocarta[1:10], type="RNASeq",  norm.method="TMM")
 PWEA(hnrnp.cnts, group, biocarta[1:10],  type="RNASeq", test="limma", nperm=100)
 TBS(hnrnp.cnts, group, biocarta[1:10], type="RNASeq", logFC.th=-1, nperm=100, test="limma")
 }
}
}
