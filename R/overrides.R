#spiaAttrs<-function(){
#out<-c("activation", "compound", "binding/association", "expression", "inhibition", "activation_phosphorylation", "phosphorylation", "inhibition_phosphorylation", "inhibition_dephosphorylation",   
# "dissociation", "dephosphorylation", "activation_dephosphorylation",   
# "state change", "activation_indirect effect", "inhibition_ubiquination",        
# "ubiquination", "expression_indirect effect", "inhibition_indirect effect",     
# "repression", "dissociation_phosphorylation", "indirect effect_phosphorylation",
# "activation_binding/association", "indirect effect", "activation_compound",            
# "activation_ubiquination")
# return(out)
#}

prepareSPIA2<-function (db, pathwaySetName, both.directions, print.names = FALSE, convert=TRUE, IDs="entrez", edgeAttrs=defaultEdgeAttrs, save.file=TRUE, reduce=TRUE)
{
    path.info <-Filter(Negate(is.null), lapply(db, function(p) {
        if (print.names) cat(p@title, "\n")
        if (convert) p <- convertIdentifiers(p, IDs)
        es <- graphite::edges(p)
        if (reduce) if (NROW(es) == 0 || NROW(unique(es[, 1:2])) < 5) return(NULL)
        ns <- graphite::nodes(p)
        spiaEdges<-edgeAttrs[[1]]
        if (!all(levels(es[,4]) %in% spiaEdges[,1] )) stop("Unexpected edge type", levels(es[,4])[!levels(es[,4]) %in% spiaEdges[,1]], "Please see '?defaultEdgeAttrs' for explanation\n")       
        es <- merge(es, spiaEdges, all.x = TRUE)[c("src", "dest", "direction", "spiaType")]
        l <- sapply(edgeAttrs[[2]][,1], simplify = FALSE, USE.NAMES = TRUE,
            function(edgeType) {
                est <- es[es[, 4] == edgeType, , drop = FALSE]
                if (both.directions) gnl <- buildGraphNEL(ns, est, TRUE) else
                 gnl<-buildGraphNEL(ns, est, FALSE)
                t(as(gnl, "matrix"))
            })
        l$title <- p@title
        l$nodes <- ns
        l$NumberOfReactions <- 0
        
        return(l)
    }))
if (save.file) save(path.info, file = paste(pathwaySetName, "SPIA.RData", sep = ""))
return(path.info)
}

runSPIA<-function(pathways, exprs, group, gene.stat="logFC", both.directions, convert=TRUE, IDs, logFC.th=2, p.val.th=0.05, edgeAttrs=defaultEdgeAttrs, ...) 
{

path.info<-prepareSPIA2(pathways,"pathways",both.directions=both.directions, convert=convert, IDs=IDs, edgeAttrs=edgeAttrs, both.directions) 


#if (!require(limma)) stop("Please, install limma package")
group<-factor(group)
if (!nlevels(group)==2) stop("Group vector has not two levels")
cat(levels(group)[1],"denoted as 0 \n", levels(group)[2],"denoted as 1\n", "Contrasts: ", levels(group)[2], "-", levels(group)[1],"\n"    )

deg.table<-testlimma(exprs, group)

de<-deg.table[abs(deg.table$logFC)>logFC.th & deg.table$P.Value<p.val.th, "logFC"]
names(de)<-deg.table[abs(deg.table$logFC)>logFC.th & deg.table$P.Value<p.val.th, "ID"]
cat("Found", length(de), "differentially expressed genes\n")
all<-deg.table$ID    

    res<-graphite.SPIA(de, all, "pathways", beta=setNames(edgeAttrs[[2]]$beta, edgeAttrs[[2]][,1]))

    names(res)[7:9]<-c("p","pFdr","pFWER")
    rownames(res)<-res[,1]
    res<-res[,-1]
    if (gene.stat=="logFC") degt<-deg.table$logFC
    if (gene.stat=="stats") degt<-deg.table$t
    names(degt)<-deg.table$ID
    
    #path.info<-prepareSPIA2(pathways,"pathways", IDs=IDs, reduce=F) 
    tsig<-spiaPF(de,all,"pathways")
    
    out<-list(res=res, topo.sig=tsig, degtest=degt)
    class(out)<-c(class(out), "topResultW","topResult")
return(out)
}





########clipper
 
runClipper<-function (pathway, expr, classes, method, both.directions, testCliques=FALSE, nperm=NULL, alphaV=0.05, ...)
{
    #if (!require(clipper))
#        stop("library clipper is missing")
#
    classes<-factor(as.numeric(factor(classes)))

    out<-sapply(pathway, function(path) {cat(path@title,"\n");
    if (both.directions) g <- buildGraphNEL(nodes(path), edges(path), TRUE) else
          g <- buildGraphNEL(nodes(path), edges(path), FALSE)
    
    
    step1<-try(pathQ(expr, classes, g, nperm=nperm, alphaV=alphaV, b=nperm, permute=TRUE),TRUE)
    step2<-try(easyClip(expr, classes, g, method = method, ...), TRUE)
    if (!inherits(step2, "try-error")) step2<-step2[which.max(step2[,"maxScore"]),c(5,7,8,11,12)] 
    if (is.list(step1) & is.data.frame(step2)) c(unlist(step1),step2) else 
    if (is.list(step1)) c(unlist(step1), rep(NA, 5)) else rep(NA,7)
        })      
      
    out<-as.data.frame(t(out))

    for (i in 1:5)     suppressWarnings(out[,i]<-as.numeric(as.character(out[,i])) )
    for (i in 6:7)     suppressWarnings(out[,i]<-as.character(out[,i])   )
    
    colnames(out)<-c("alphaVar","alphaMean", "maxScore","activation","impact","involvedGenes", "pathGenes")

    
    if (testCliques) {
    cat("Testing cliques...\n")
    cat("------------------\n")
    cli<-lapply(pathway, function(path) {cat(path@title,"\n");
      if (both.directions) g <- buildGraphNEL(nodes(path), edges(path), TRUE) else
          g <- buildGraphNEL(nodes(path), edges(path), FALSE)

     if (method=="mean") cliqueMeanTest(expr, classes, g, nperm, alphaV=alphaV, b=nperm, root=NULL, permute=TRUE) else
                       cliqueVarianceTest(expr, classes, g, nperm, alphaV=alphaV, b=nperm, root=NULL, permute=TRUE)
      })


    }else cli<-NULL
    out<-list(res=out, topo.sig=cli, degtest=NULL)
    class(out)<-c(class(out),  "topResultC","topResult")
    return(out)
}

#TopologyGSA
runTopologyGSADirected<-function(pathway, test, exp1, exp2, alpha, both.directions, ...)
{
   # if (!require(topologyGSA))
#        stop("library topologyGSA is missing")
#        checkPkgVersion<-function (name, min_version)
#{
#    version <- package_version(installed.packages()[name, "Version"])
#    if (version < package_version(min_version))
#        stop("the installed ", name, " version is too old (need at least ",
#            min_version, ")")
#}
#    checkPkgVersion("topologyGSA", "1.0")
   if (insufficientCommonGenes(pathway, colnames(exp1))) 
        return(NULL)
   test<-switch(test, var = pathway.var.test, mean = pathway.mean.test, 
        stop("invalid test type: ", name))

       if (both.directions) g <- buildGraphNEL(nodes(pathway), edges(pathway), TRUE) else
          g <- buildGraphNEL(nodes(pathway), edges(pathway), FALSE)
    test(exp1, exp2, g, alpha, ...)
}




