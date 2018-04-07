#' Pathway regulation score (PRS) 
#' 
#' This function implements the PRS method to analyze pathway enrichment of gene
#' expression data. For PRS, a gene weight correspond to the number of 
#' downstream differentially expressed genes. 
#'
#' 
#' @aliases prsWeights
#' @param de A named numeric vector containing log2 fold-changes of the 
#' differentially expressed genes. Recommended names are Entrez gene IDs. 
#' @param all A character vector with the gene IDs in the reference set. If the
#' data was obtained from a gene expression experiment, this set will
#' contain all genes measured in the  experiment. This vector should contain 
#' *all* names of the \code{de} argument.
#' @param pwy A \code{linkS4class{Pathway}} for which the weights should be
#' computed. 
#' @param pwys A \code{linkS4class{PathwayList}} containing the pathways that 
#' should be analyzed for enrichment. 
#' @param nperm Integer. Number of permutations.
#' @return A \code{data.frame} with normalized score and p-value for each 
#' pathway analyzed.
#' @author Ivana Ihnatova
#' @seealso \code{\link{pathways}} 
#' @references Ibrahim et al. (2012) A topology-based score for pathway 
#' enrichment. J Comput Biol, 19(5):563-73.
#' @examples
#'
#' # pathways
#' library(graphite)
#' pwys <- pathways("hsapiens","kegg")[1:10]
#' 
#' # expression data
#' all <- nodes(pwys[[1]])
#' nds <- sample(all, 30)
#' de <- setNames(rnorm(30), nds)
#' 
#' # executing PRS
#' prsWeights(pwys[[1]], de, all)
#' prs(de, all, pwys, nperm=100) 
#' 
#' @export
prs <- function(de, all, pwys, nperm=1000) 
{
    # remove ID type
    if(all(grepl(":", all)))
    {
        .spl2 <- function(s) unlist(strsplit(s, ":"))[2]
        all <- vapply(all, .spl2, character(1), USE.NAMES=FALSE)
        names(de) <- vapply(names(de), .spl2, character(1), USE.NAMES=FALSE)
    }    

    de <- abs(de) #yes?    
    perms <- .preparePerms(de, all, nperm, method = "PRS")
    pwys <- .preparePathways(pwys, method="PRS", genes=all)
    res <- vapply(pwys, 
                    function(p) .PRSSingle(p, de, all, perms),
                    numeric(2))
    
    res <- data.frame(t(res))
    sorting.df <- res[,c(2,1)]
    sorting.df[,2] <- -abs(sorting.df[,2]) 
    res <- res[do.call(order, sorting.df),]
    return(res)
}

#' @rdname prs
#' @export
prsWeights <- function(pwy, de, all) 
{
    # remove ID type
    if(all(grepl(":", all)))
    {
        .spl2 <- function(s) unlist(strsplit(s, ":"))[2]
        all <- vapply(all, .spl2, character(1), USE.NAMES=FALSE)
        names(de) <- vapply(names(de), .spl2, character(1), USE.NAMES=FALSE)
    }

    if(is(pwy, "Pathway")) 
        pwy <- .preparePathways(list(pwy), method="PRS", genes=all)[[1]]

    g <- intersect(all, rownames(pwy))
    glen <- length(g)

    if(!glen) stop("Pathway does not contain any measured genes")
    set <- pwy[g, g]
    
    wei <- setNames(rep(0, glen), g)
    ind <- g %in% names(de)
    if(sum(ind))
    { 
        weight <- 1 + downstreamCpp(set, g, g[ind])
        wei[names(weight)] <- weight
    }
    return(wei)
}

.PRSSingle <- function(pwy, de, all, perms) 
{
    weight <- prsWeights(pwy, de, all)
    
    g <- intersect(all, rownames(pwy))
    glen <- length(g)
    set <- pwy[g, g]
    indde <- g %in% names(de)

    nf <- sum(indde) / glen
    expr <- rep(1, glen)
    expr[indde] <- de[g[indde]] 
    
    obs <- sum(expr * weight) * nf
    
    # random
    permsub <- perms[g, , drop = FALSE]
    
    weight.rn <- apply(permsub, 2, function(x)
    {
        weight <- setNames(rep(0, glen), g)
        if (glen && length(g[x != 0]))
        { 
            lind <- g[x != 0]
            weight[lind] <- 1 + downstreamCpp(set, g, lind)        
            weight[x == 0] <- 0
        }
        return(weight)
    })
    
    nf.rn <- colSums(permsub != 0)
    rand <- colSums(weight.rn * permsub) * (nf.rn / glen)
    
    # normalization
    obs <- (obs - mean(rand)) / sd(rand)
    rand <- (rand - mean(rand)) / sd(rand)
    p.value <- (sum(rand >= obs) + 1) / (length(rand) + 1)
    res <- c(nPRS = obs, p.value = p.value)
    return(res)
}


