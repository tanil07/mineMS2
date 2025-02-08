###The ms2LibSplit is similar to the ms2Lib object, except that it is made to compute patterns on components.


##Checking the correctness of the found components.
checkComponent <- function(m2l,comp){
	if(is.numeric(comp)|is.integer(comp)){
		valid <- all(comp<length(mm2Spectra(m2l)))
		return(list(comp=unique(comp),valid=valid))
	}

	if(all(startsWith(comp,"S"))){
		return(list(comp=unique(as.numeric(str_sub(comp,2))),valid=valid))
	}

	###Checking if they corresponds to IDs.
	idm <- match(comp,mm2Ids(m2l))
	if(any(is.na(idm))){
		return(list(comp=comp,valid=FALSE))
	}else{
		return(list(comp=unique(idm),valid=TRUE))
	}
	return(list(comp=comp,valid=FALSE))
}




#' Set the components on an ms2LibSplit object.
#'
#' Set the compoentns field of an ms2LibSplit object, whilme checking for their correctness.
#'
#' @param m2l An ms2LibSplit object.
#' @param components A list containing the components. Each components needs to be one of the following
#' \itemize{
#' \item A vector corresponding to furnished spectra index.
#' \item A vector of ids corresponding to the ids of the spectra.
#' \item A vecotr of spectra ids starting with S.
#' }
#'
#' @return The filled ms2Lib object.
#' @export
#'
#' @examples
#' print("examples to be put here")
setMethod("setComponents","ms2LibSplit",function(m2l,components){

	###Checking the validity of the furnished compoents.
	components <- sapply(components,checkComponent)

	valids <- sapply(components,function(x) x$valid)
	
	if(any(!valids)){
		stop("Invalids ids in component(s): ",paste(components[which(!valids)],sep=", "))
	}

	m2l@components <- sapply(components,function(x) x$comp)
	m2l
})


setMethod("mineClosedSubgraphsByComponent","ms2LibSplit",function(m2l, count = 2, sizeMin = 2, precursor = FALSE){

	if(length(m2l@components)==0) stop("Components should be set using the setComponents methods before the mining for an ms2LibSplit object.")

	if(count<2){
		warning("'count' parameter too low; value set to 2.")
		count <- 2
	}


	###Get the data.frame corresponding to the sizes.
	processing <- sapply(mm2Dags(m2l),function(x){
		ecount(x)>1
	})

	if(nrow(mm2EdgesLabels(m2l))==0){
		stop("No labels constructed, use the discretizeMzDifferences method first.")
	}

	if(sizeMin==1&nrow(mm2EdgesLabels(m2l))>600){
		###Wide variety of mass losses.
		warning("Signle edges graphs allowed, risk of computational overhead.")
	}

	###All the dags are converted before processing.
	df_edges <- sapply(mm2Dags(m2l),fromIgraphToDf_edges,simplify = FALSE)
	df_vertices <-sapply(mm2Dags(m2l),fromIgraphToDf_vertices,simplify = FALSE)

	kTree <- NULL
	if(length(mm2EdgesLabels(m2l))<600){
		kTree <- 2
	}else{
		kTree <- 1
	}

	allPatterns <- vector(mode="list",length=length(m2l@components))

	tempidx <- numeric(length(m2l@components)+1)
	tempidx[1] <- 0

	###Now we start the detection component by component.
	for(icomp in seq_along(m2l@components)){
		comp <- m2l@components[[icomp]]
		###Mining the patterns.
		resRcpp <- mineClosedDags(df_vertices[comp],df_edges[comp],
								  processing,count,kTree,sizeMin,precursor)

		###We build the fragmentation pattern for each components.
		res <- sapply(resRcpp$patterns[comp],fragPattern,USE.NAMES =  FALSE)

		res <- sapply(res, relabelOccurrences,newlabs=comp)

		allPatterns[[i]] <- res

		tempidx[i] <- length(res)
	}



	# mm2LatticeIdxItems(m2l) <- resRcpp$items

	###Construction via fragPatternc constructor.
	mm2Patterns(m2l) <- do.call('c',allPatterns)

	###Initializing the names of the patterns.
	for(i in 1:length(m2l@patterns)) mm2Name(m2l@patterns[[i]]) <- paste("P",i,sep="")


	###Buillding of the lattice
	# mm2Lattice(m2l) <- graph_from_data_frame(resRcpp$edges,directed=TRUE,resRcpp$nodes)

	###We add a label filed which give all the values of this.

	message("Processing finished, ",length(mm2Patterns(m2l))," patterns mined.")
	m2l
})

#' Show an ms2LibSplit object.
#'
#' Get the string representation of an ms2LibSplit object
#'
#' @param object An m2Lib object to be shown.
#'
#' @return None.
#' @export
setMethod("show","ms2LibSplit",function(object){
	cat("An ms2Lib object containing",length(object),"spectra.\n")
	cat("It has",nrow(mm2EdgesLabels(object)),"edges labels.\n")
	cat("The available supplementary informations are:",colnames(mm2SpectraInfos(object)),"\n")
	cat("It contains:",length(mm2Patterns(object)),"patterns\n")
	if(length(mm2ReducedPatterns(object))!=0) cat("It has been reduced to",
											   length(mm2ReducedPatterns(object)),"patterns")
})
