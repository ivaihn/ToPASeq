spiapf<-function (de = NULL, all = NULL, organism = "hsa", data.dir = NULL,
    pathids = NULL, nB = 2000, verbose = TRUE,
    beta = NULL, combine = "fisher")
{
    if (is.null(de) | is.null(all)) {
        stop("de and all arguments can not be NULL!")
    }
    rel <- c("activation", "compound", "binding/association",
        "expression", "inhibition", "activation_phosphorylation",
        "phosphorylation", "inhibition_phosphorylation", "inhibition_dephosphorylation",
        "dissociation", "dephosphorylation", "activation_dephosphorylation",
        "state change", "activation_indirect effect", "inhibition_ubiquination",
        "ubiquination", "expression_indirect effect", "inhibition_indirect effect",
        "repression", "dissociation_phosphorylation", "indirect effect_phosphorylation",
        "activation_binding/association", "indirect effect",
        "activation_compound", "activation_ubiquination")
    if (is.null(beta)) {
        beta = c(1, 0, 0, 1, -1, 1, 0, -1, -1, 0, 0, 1, 0, 1,
            -1, 0, 1, -1, -1, 0, 0, 1, 0, 1, 1)
        names(beta) <- rel
    } else {
        if (!all(names(beta) %in% rel) | length(names(beta)) !=
            length(rel)) {
            stop(paste("beta must be a numeric vector of length",
                length(rel), "with the following names:", "\n",
                paste(rel, collapse = ",")))
        }
    }
    .myDataEnv <- new.env(parent = emptyenv())
    datload <- paste(organism, "SPIA", sep = "")
    load(file = paste(getwd(), paste(datload, ".RData",
                sep = ""), sep = "/"), envir = .myDataEnv)
        
    datpT = .myDataEnv[["path.info"]]
    if (!is.null(pathids)) {
        if (all(pathids %in% names(datpT))) {
            datpT = datpT[pathids]
        }   else {
            stop(paste("pathids must be a subset of these pathway ids: ",
                paste(names(datpT), collapse = " "), sep = " "))
        }
    }
    datp <- list()
    path.names <- NULL
    hasR <- NULL
    for (jj in 1:length(datpT)) {
        sizem <- dim(datpT[[jj]]$activation)[1]
        s <- 0
        con <- 0
        for (bb in 1:length(rel)) {
            con = con + datpT[[jj]][[rel[bb]]] * abs(sign(beta[rel[bb]]))
            s = s + datpT[[jj]][[rel[bb]]] * beta[rel[bb]]
        }
        z = matrix(rep(apply(con, 2, sum), dim(con)[1]), dim(con)[1],
            dim(con)[1], byrow = TRUE)
        z[z == 0] <- 1
        datp[[jj]] <- s/z
        path.names <- c(path.names, datpT[[jj]]$title)
        hasR <- c(hasR, datpT[[jj]]$NumberOfReactions >= 1)
    }
    names(datp) <- names(datpT)
    names(path.names) <- names(datpT)
    tor <- lapply(datp, function(d) {
        sum(abs(d))
    }) == 0 | hasR | is.na(path.names)
    #datp <- datp[!tor]
    #path.names <- path.names[!tor]
    IDsNotP <- names(de)[!names(de) %in% all]
    if (length(IDsNotP)/length(de) > 0.01) {
        stop("More than 1% of your de genes have IDs are not present in the reference array!. Are you sure you use the right reference array?")
    }
    if (!length(IDsNotP) == 0) {
        cat("The following IDs are missing from all vector...:\n")
        cat(paste(IDsNotP, collapse = ","))
        cat("\nThey were added to your universe...")
        all <- c(all, IDsNotP)
    }
    if (length(intersect(names(de), all)) != length(de)) {
        stop("de must be a vector of log2 fold changes. The names of de should be included in the refference array!")
    }
    ph <- pb <- pcomb <- NULL
    set.seed(1)
    res<-list()
 
    for (i in 1:length(names(datp))) {
        if (tor[i]) {Acc<-c();logFC<-c()} else {
        path <- names(datp)[i]
        #cat(path,"\n")
        M <- datp[[path]]
        diag(M) <- diag(M) - 1
        X <- de[rownames(M)]
        noMy <- sum(!is.na(X))
        okg <- intersect(rownames(M), all)
        ok <- rownames(M) %in% all
        Acc<-c()
        logFC<-c()
        if ((noMy) > 0 & (abs(det(M)) > 1e-07)) {
            X[is.na(X)] <- 0
            pfs <- solve(M, -X)
            Acc=X
            logFC=pfs-X
       }}
       res[[i]]<-t(data.frame(Acc, logFC)  )
    }
    names(res)<-names(datp)
    res<-lapply(res, function(x) if(ncol(x)>0) x else NULL)
    res
}

spiaPF<-function (de, all, pathwaySetName, ...) 
{
    optArgs <- list(...)
    optArgs$organism <- pathwaySetName
    spiaEnv <- environment(spia)
    fakeEnv <- new.env(parent = spiaEnv)
    fakeFile<-function(name, package=NULL, ...){
    if (package == "SPIA" && name == "extdata") 
        getwd() else fakeFile<-base::system.file(name, package = package, ...)
        }
        
    assign("system.file", fakeFile, fakeEnv)
    tryCatch({
        environment(spia) <- fakeEnv
        do.call(spiapf, c(list(de, all), optArgs))
    }, finally = {
        environment(spia) <- spiaEnv
    })
}
