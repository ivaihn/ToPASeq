setMethod("convertIdentifiers", c("pathway","character"),
  function(x, to) {

    db <- loadDb(x@species)
    to <- destIdentifier(to, db)

    if (x@ident != to) {

      cn <- colnames(x@edges)

      tblsrc <- convert(db, x@edges, "src", x@ident, to)
      tbldest <- convert(db, x@edges, "dest", x@ident, to)
      tbl<-c(tblsrc, tbldest)
      tbl<-tbl[unique(names(tbl))]
      convertIdentifiersByVector(x, tbl)
      x@ident <- to

      if (nrow(x@edges) == 0)
        warning("the conversion lost all edges of pathway \"", x@title, "\"")
    }

    return(x)
  })


loadDb <- function(species) {
  db <- selectDb(species)

  tryCatch(get(db),
    error=function(e) {
      if (!suppressPackageStartupMessages(require(db, character.only=TRUE)))
        stop("package \"", db, "\" is missing", call.=FALSE)
      get(db)
    })
}

destIdentifier <- function(to, db) {

  if (!(to %in% columns(db)))
    stop(to, " is not supported with this species. 
 Available types are ",paste(columns(db), collapse=" "),
         call.=FALSE)

  return(to)
}

selectDb <- function(species) {

  l <- list(athaliana="org.At.tair.db",
            btaurus="org.Bt.eg.db",
            celegans="org.Ce.eg.db",
            cfamiliaris="org.Cf.eg.db",
            dmelanogaster="org.Dm.eg.db",
            drerio="org.Dr.eg.db",
            ecoli="org.EcK12.eg.db",
            ggallus="org.Gg.eg.db",
            hsapiens="org.Hs.eg.db",
            mmusculus="org.Mm.eg.db",
            rnorvegicus="org.Rn.eg.db",
            scerevisiae="org.Sc.sgd.db",
            sscrofa="org.Ss.eg.db",
            xlaevis="org.Xl.eg.db")

  n <- l[[species]]
  if (is.null(n))
    stop("no such species")

  return(n)
}

convert <- function(db, edges, colname, from, to) {

  keys <- unique(edges[[colname]])
  tbl <- lookupKeys(db, keys, from, to)
  tbl <- setNames(tbl[,2], tbl[,1])
  return(tbl)
}

lookupKeys <- function(db, keys, from, to)
  tryCatch(
    na.omit(suppressWarnings(select(db, keys, to, from))),
    error=function(e) NULL)



