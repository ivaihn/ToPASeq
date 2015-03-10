#conv.table = named vector of new identifiers, names refer to the original ones
#id.type = type of the new identifiers  e.g. "TAIR"

convertIdentifiersByVector<-function (pathway, conv.table, id.type="unknown") 
{
    ns <- nodes(pathway)
    es <- edges(pathway)
    if (!all(ns %in% names(conv.table))) 
      warning(paste("These pathway nodes are missing in the 'conv.table':", paste(ns[! ns %in% names(conv.table)], collapse=", "),". The original identifiers were kept."))


    miss <- ns[! ns %in% names(conv.table)]
    conv.table<-c(conv.table, setNames(miss,miss))
    ns <- unname(conv.table[nodes(pathway)])
    es[,1]<- unname(conv.table[es[,1]])
    es[,2]<- unname(conv.table[es[,2]])

    pathway@nodes <- ns
    pathway@edges <- es
    pathway@ident <- id.type
    return(pathway)
}

