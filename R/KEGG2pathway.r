KEGG2Pathway<-function(file, expandGenes=TRUE, expandCom=TRUE, nongene=c("keep","propagate", "discard"), 
 ident="KEGGnative", database="KEGG",
 species=NULL){
kegg.path<-parseKGML(file) 

path<-kegg.path@pathwayInfo
path.df<-data.frame(name=path@name, org=path@org, number=path@number,
  title=path@title, image=path@image, link=path@link)
nodes<-kegg.path@nodes
nodes.df<-sapply(nodes, function(x) {
  if (x@type=="group") na<-x@component else na<-x@name
  data.frame(entryID=as.character(x@entryID),
  name=paste(as.character(na),collapse=", "), type=as.character(x@type), link=as.character(x@link),
  reaction=as.character(x@reaction), map=as.character(x@map),
  graph.name=as.character(x@graphics@name), graph.x=x@graphics@x,
  graph.y=x@graphics@y, graph.type=as.character(x@graphics@type), 
  graph.width=x@graphics@width, graph.height=x@graphics@height, 
  graph.fgcolor=as.character(x@graphics@fgcolor),
  graph.bgcolor=as.character(x@graphics@bgcolor),stringsAsFactors=FALSE )}
  )

colnames(nodes.df)<-NULL
interactions<-kegg.path@edges

if (length(interactions)> 0) { 

interactions.df<-sapply(interactions, function(x) {
 if (length(x@subtype)>0) {
  subtype.name<-paste(sapply(x@subtype, function(y) y@name), collapse=", ")
  subtype.value<-paste(sapply(x@subtype, function(y) y@value), collapse=", ")  
  } else {
  subtype.name<-subtype.value<-NA
  }
 c( 
 entry1ID=x@entry1ID, entry2ID=x@entry2ID, type=x@type, 
 subtype.name=subtype.name, subtype.value=subtype.value)
})
colnames(interactions.df)<-NULL


reactions<-kegg.path@reactions

if (is.null(species)) species<-substring(file,1,regexpr("[:digits:]", "./hsa01234.xml")-1)

nod<-t(nodes.df)[,1:3]
if (expandGenes) {
exp.count<-sapply(nod[,2], function(x) length(strsplit(x,", ")[[1]]))
exp.names<-unlist(sapply(nod[,2], function(x) strsplit(x,", ")[[1]]))
exp.ids<-unlist(Map(function(x,y) rep(x,y),nod[,1], as.list(exp.count)))
exp.type<-unlist(Map(function(x,y) rep(x,y),nod[,3], as.list(exp.count)))
nod<-cbind(exp.ids, exp.names, exp.type)
} else {
colnames(nod)<-c("exp.ids","exp.names","exp.type")
}

if (sum(nod[,3]=="gene")==0) {
message("The pathway does not contain any genes. Returning empty pathway")
return(NULL)
}

nod[nod[,3]=="group",2]<-nod[match(nod[nod[,3]=="group",2], nod[,1]),2]


if (nongene=="discard") {nod<-nod[nod[,3]=="gene"| nod[,3]=="group",]} 


interactions.df<-t(interactions.df)


cmi<-interactions.df[interactions.df[,4] %in% c("compound", "hidden compound"),, drop=FALSE]
cmi<-rbind(cmi[,c(1,5,3,4,5)], cmi[,c(5,2,3,4,5)])
cmi[,5]<-""

interactions.df<-interactions.df[! interactions.df[,4] %in% c("compound","hidden compound"),, drop=FALSE]
interactions.df<-rbind(interactions.df, cmi)


edg.tabs<-lapply(seq_len(nrow(interactions.df)), function(i) {
      x<-interactions.df[i,]
      src<-nod[nod[,1]==x[1],2]
      des<-nod[nod[,1]==x[2],2]
      #if (union(src, des)!=intersect(src,des)) {
      out<-rbind(expand.grid(setdiff(src,des), des, unname(x[3]),unname(strsplit(x[4],", ")[[1]]), unname(x[5])),
                 expand.grid(intersect(src,des), setdiff(des, src),   unname(x[3]),unname(strsplit(x[4],", ")[[1]]), unname(x[5])))
      #} else out<-expand.grid(src, des, unname(x[3]), unname(strsplit(x[4],", ")[[1]]), unname(x[5])) 
      attr(out,"out.attrs")<-NULL
      colnames(out)<-c("src","dest", "type", "name", "value")
      rownames(out)<-NULL
      return(out)
    })
    
edg.tabs<-Reduce(rbind,edg.tabs)

if (nrow(edg.tabs)==0) {
message("Pathway contains only self loops. Returning NULL")
return(NULL)
}
levels(edg.tabs[,4])[levels(edg.tabs[,4])=="compound"]<-"indirect effect"

E<-data.frame(src=as.character(edg.tabs[,1]), dest=as.character(edg.tabs[,2]), 
 direction=factor(rep("directed",nrow(edg.tabs))), type=factor(paste("process(",edg.tabs[,4],")", sep="")), stringsAsFactors=FALSE)


if (expandCom) {
groups<-nod[nod[,3]=="group",]
comp<-tapply(groups[,2], groups[,1], function(x) x)
compEd<-lapply(comp, function(x) {
 edMat<-t(combn(x,2))
 n<-nrow(edMat)
 data.frame(src=edMat[,1], dest=edMat[,2], direction=rep("undirected",n), type=rep("binding",n ), stringsAsFactors=FALSE )
 })
 compEd<-Reduce(rbind,compEd)
 E<-rbind(E,compEd)
}

E<-E[!duplicated(E),]

gr<-new("Pathway", id=path@name, title=path@title,  edges=E, 
database=database,species=species, identifier=ident,  
timestamp=Sys.Date())


if (nongene=="propagate") {
 att<-sapply(nodes(gr), function(x) strsplit(x,":",fixed=TRUE)[[1]][1])
 
 non.genes<-nodes(gr)[att %in% c("path","ko", "ec", "rn","cpd","gl","group")] #alebo len cisla

  for (el in non.genes) gr<-eliminateNode(gr, el)
 }
} else {
message("The pathway does not contain any edges. Returning NULL")
return(NULL)
}
return(gr) 
}

KEGG2pathway<-function(file, expandGenes=TRUE, expandCom=TRUE, nongene=c("keep","propagate", "discard"), 
 ident="KEGGnative", database="KEGG",
 species=NULL){
 .Deprecated("KEGG2Pathway")
 }