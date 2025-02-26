#' @include references.R

## @name mm2Graph
## @aliases fragPattern-mm2Graph


#' Methods for a fragPattern objects 
#'
#' @rdname fragPattern-methods
#' @export
#' @examples 
#' data(m2l)
#' ## igraph object to store the vertices and edges of a pattern
#' graph <- mm2Graph(m2l["P1"])
setMethod("mm2Graph", "fragPattern",function(object){
	return(object@graph)
})


#' @rdname fragPattern-methods
#' 
#' @export
#' @examples 
#' 
#' occs <- mm2Occurrences(m2l["P1"])
setMethod("mm2Occurrences", "fragPattern",function(object){
	return(object@occurrences)
})

#' @rdname fragPattern-methods
#' 
#' @export
#' @examples 
#' 
#' root <- mm2Root(m2l["P1"])
setMethod("mm2Root", "fragPattern",function(object){
	return(object@root)
})

#' @rdname fragPattern-methods
#' 
#' @export
#' @examples 
#' 
#' name_p <- mm2Name(m2l["P1"])
setMethod("mm2Name", "fragPattern",function(object){
	return(object@name)
})

#' @rdname fragPattern-methods
#' 
#' @export
#' @examples 
#' 
#' can_form <- mm2CanonicalForm(m2l["P1"])
setMethod("mm2CanonicalForm", "fragPattern",function(object){
  return(object@canonicalForm)
})

setMethod("mm2Graph<-", "fragPattern", function(pat,value){
	pat@graph <- value
	pat
})

setMethod("mm2Occurrences<-", "fragPattern", function(pat,value){
	pat@occurrences <- value
	pat
})

setMethod("mm2Root<-","fragPattern", function(pat,value){
	pat@root <- value
	pat
})

setMethod("mm2Name<-", "fragPattern",function(pat,value){
	pat@name <- value
	pat
})

setMethod("mm2CanonicalForm<-", "fragPattern",function(pat,value){
  pat@canonicalForm <- value
  pat
})


#' Show information about a fragPattern object
#'
#' @param object A fragPattern object
#' @export
#'
#' @rdname fragPattern-methods
#' 
#' @return None
#' @examples
#' 
#' show(m2l["P12"])
setMethod("show","fragPattern",function(object){
	cat("A fragPattern object with",vcount(mm2Graph(object))-1,"losses occuring ",
		nrow(mm2Occurrences(object)),"times.\n")
})


#' fragPattern constructor
#'
#' This is an internal constructor which is used to pass from the rcpp
#' output to a fragPattern object including the occurrences. It should not be
#' called directly by the user.
#'
#' @param rlist rlist is a list containg 2 fields,
#' \code{edges}: An edge list which includes the "from","to" and "lab" fields, which be processed
#' by the graph_from_data_frame function.
#' \code{occurrences}: A matrix with 2 columns, the first giving the id in which the patterns are found
#' and the second giving the nodes.
#'
#' @return
#' The constructed fragPattern object.
fragPattern <- function(rlist){
	fp <- new("fragPattern")
	mm2Graph(fp) <- graph_from_data_frame(rlist$edges)
	mm2Occurrences(fp) <- rlist$occurrences

	##The root is always 0
	mm2Root(fp) <- as.integer(0)
	fp
}


###Make an ID fo the graph based on edges.
makeId <- function(pat){
  elabs <- edge_attr(pat@graph,'lab',incident_edges(pat@graph,1)[[1]])
  elabs <- paste(sort(elabs),collapse="/")
  elabs
 }

###Building the canonical form of a pattern as an itemsep
canonicalForm <- function(pat){
  root_labels <- edge_attr(pat@graph,"lab",incident_edges(pat@graph,1)[[1]])
  pat@canonicalForm <- paste(root_labels,collapse = "|")
  return(pat)
}


#' Calculate the coverage of a specific pattern
#' 
#' The S4 method for the fragPattern object should never be called by the user. Call the ms2Lib method instead.
#'
#' @param x The pattern
#' @param mzloss The table of m/z differences
#' @param mgs the m/z graphs (the DAGs)
#' @export 
#' 
#' @rdname calculateCoverage
#'
#' @return The m2l object with all the coverage calculated.
setMethod("calculateCoverage","fragPattern",function(x,mzloss,mgs){
  occs <- x@occurrences
  coverages <- rep(NA_real_,nrow(occs))
  for(j in 1:nrow(occs)){
    gid <- occs[j,"gid"]
    cg <-mgs[[gid]]
    mapv <- get_mapping(cg ,x@graph, mzloss,root=occs[j,"idx"])
    
    ###Plotting of the spectra
    intv <- vertex_attr(cg, "rel_int")
    mzs <- vertex_attr(cg, "mz")
    ids <- V(cg)
    
    ###Peaks are split between matched and non matched.
    matched_peaks_idx <- match(mapv[2,], ids)
    coverages[j] <- sum(intv[matched_peaks_idx])/sum(intv)
  }
  onames <- colnames(occs)
  if(!COVERAGE_NAME %in% colnames(occs)){
    occs <- cbind(occs,coverages)
    colnames(occs) <- c(onames,COVERAGE_NAME)
  }else{
    occs[,COVERAGE_NAME] <- coverages
  }
  x@occurrences <- occs
  x
})

setMethod("hasCoverage","fragPattern",function(x){
  COVERAGE_NAME %in% colnames(x@occurrences)
})


relabelOccurrences <- function(pat,newlabs){
	tempOccs <- mm2Occurrences(pat)
	tempOccs[,1] <- newlabs[tempOccs[,1]]
	mm2Occurrences(pat) <- tempOccs
	pat
}

#' Spectra information for a specific pattern
#' 
#' @param m2l a ms2Lib object
#' @param id_p a fragPattern id (such as "P1")
#' 
#' @return a data.frame containing information about the spectra containing the pattern
#' @export
infoPatterns <- function(m2l, id_p)
{
	occs <- mm2Occurrences(m2l[id_p])[,'gid']
	s_occs <- sapply(occs, function(x){return(paste0("S", x))})
	infos <- getInfo(m2l, "S")

	df <- infos[rownames(infos) %in% s_occs, names(infos) != "file"]
  return(df)
}



