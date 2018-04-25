#' Plot patterns and occurences of an ms2Lib object
#'
#' Plot patterns graph and occurences of an ms2Lib object.
#'
#' @param m2l An ms2Lib object
#' @param ids THe ids to be plotted or NULL if all the ids are supposed to be plot
#' @param occurences include the occurences
#' @param full Shall the full patterns be plotted (it can take some times)
#' @param byPage Maximum number of occurences by page.
#'
#' @export
#'
#' @examples
plotPatterns <- function(m2l,ids=NULL,occurences=TRUE,full=FALSE,byPage=9,...){


	if(is.null(ids)){
		if(length(mm2ReducedPatterns(m2l))==0){
			if(full){
				ids <- vrange(m2l,"P",reduced=FALSE)
			}else{
				stop("No reduced pattern, to plot the full patterns set the 'full' parameter to TRUE")

			}
		}else{
			ids <- vrange(m2l,"P",reduced=TRUE)
		}
	}else{
		tempids <- parseId(ids)
		seq_type <- sapply(tempids,"[[",i=1)
		if(any(seq_type != "P")){
			stop("Invalid ids furnished ",ids[which(seq_type!="P")]," only patterns ids may be furnished")
		}
	}

	for(id in ids){
		plot(m2l,id)
		plotOccurences(m2l,id)
	}
}

#' Plot the lattice of an ms2Lib object
#'
#' @param m2l The ms2Lib object
#' @param reduced Shall the reduced lattice or the full lattice be plot.
#'
#' @return
#' @export
#'
#' @examples
plotLattice <- function(m2l,reduced=TRUE){
	if(vcount(mm2Lattice(m2l))==0) stop("Patterns should have been mined to plot the lattice.")
	if(reduced){
		if(vcount(mm2ReducedLattice(m2l))==0) stop("Reduced lattice should have been constructed.")
	}else{

		if(vcount(mm2Lattice(m2l))==0) stop("Lattice should have been constructed.")
		stop("Not implemented yet.")

	}


}
