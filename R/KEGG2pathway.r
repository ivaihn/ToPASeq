KEGG2pathway<-function(file, expandGenes=TRUE, expandCom=TRUE, nongene=c("keep","propagate", "discard"), 
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

if (is.null(species)) species<-substring(file,1,3)

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
return(new("pathway", title=path@title, nodes=character() , 
edges=data.frame(src=character(), dest=character(), direction=factor(),
 type=factor(), stringsAsFactors=FALSE), 
ident=ident, database=database, species=species, timestamp=Sys.Date()  ))
}

nod[nod[,3]=="group",2]<-nod[match(nod[nod[,3]=="group",2], nod[,1]),2]


if (nongene=="discard") {nod<-nod[nod[,3]=="gene"| nod[,3]=="group",]} 


interactions.df<-t(interactions.df)
edg.tabs<-lapply(seq_len(nrow(interactions.df)), function(i) {
x<-interactions.df[i,]
ns<-length(nod[nod[,1]==x[1],2])
nd<-length(nod[nod[,1]==x[2],2])
ni<-length(strsplit(x[4],", ")[[1]])
src<-rep(nod[nod[,1]==x[1],2],each=nd*ni) 
des<-rep(nod[nod[,1]==x[2],2],times=ns*ni) 
ty<-rep(x[3], ns*nd)
na<-rep(rep(strsplit(x[4],", ")[[1]], each=nd),times=ns)
va<-rep(x[5], ns*nd)
out<-data.frame(src=src, dest=des, type=ty, name=na, value=va, row.names = NULL)

return(out)
})

edg.tabs<-Reduce(rbind,edg.tabs)

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


p<-new("pathway", title=path@title, nodes=unique(nod[,2]) , edges=E, 
ident=ident, database=database, species=species, timestamp=Sys.Date())
return(p)


if (nongene=="propagate") {
 att<-sapply(nodes(p), function(x) strsplit(x,":",fixed=TRUE)[[1]][1])
 
 non.genes<-nodes(p)[att %in% c("path","ko", "ec", "rn","cpd","gl","group")] #alebo len cisla
 p<-eliminateNode(p, non.genes)
 }
 
}