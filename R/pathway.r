setClass("pathway",
         representation(title="vector",
                        nodes="vector",
                        edges="data.frame",
                        ident="vector",
                        database="vector",
                        species="vector",
                        timestamp="Date"))

setMethod("show", signature(object="pathway"),
          function(object) {
            cat(paste('"', object@title, '" pathway from ', object@database, "\n",
                      "Number of nodes     = ", length(object@nodes), "\n",
                      "Number of edges     = ", NROW(object@edges), "\n",
                      "Type of identifiers = ", object@ident, "\n",
                      "Retrieved on        = ", object@timestamp, "\n", 
                      sep=""))
          })

setMethod("nodes", signature(object="pathway"), function(object) object@nodes)

setMethod("edges", signature(object="pathway"), function(object) object@edges)

setAs("Pathway","pathway", function(from) {
E<-from@edges
E[,4]<-factor(E[,4])
out<-new("pathway", title =  from@title, nodes = as.character(unname(unlist(from@edges[,1:2]))), edges = E, 
 ident = from@identifier, database =  from@database, species = from@species, timestamp = from@timestamp)
return(out)
})

setAs("pathway","graphNEL", function(from) {
out<-buildGraphNEL(nodes(from), edges(from), TRUE)
})

