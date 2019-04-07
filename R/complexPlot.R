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
plotPatterns <- function(m2l,ids=NULL,components = NULL,occurences=TRUE,full=FALSE,byPage=9,titles=NULL,...){


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
		tempids <- sapply(ids,parseId,m2l=m2l,simplify = FALSE)
		seq_type <- sapply(tempids,"[[",i="type")
		if(any(seq_type != "patterns")){
			stop("Invalid ids furnished ",ids[which(seq_type!="P")]," only patterns ids may be furnished")
		}
	}

	if(!is.null(titles) & (length(titles)!= length(ids))){
		stop("titles argument should be a chracter of the same size than ids.")
	}

	if(!is.null(components) & (length(components)!= length(ids))){
		stop("components argument should be a list of the same size than ids.")
	}

	for(pid in seq_along(ids)){
		id <- ids[pid]
		if(!is.null(titles)){
			if(!is.null(components)){
				layout(matrix(c(1,2,2,2),byrow=TRUE,ncol=4))
				plot(0,xaxt="n",yaxt="n",xlab="",ylab="",ylim=c(0,1),xlim=c(0,1),type="n",bty="n")
				text(x = rep(0,length(components[[pid]])),y=seq(0.1,0.9,length=length(components[[pid]])),
					 labels = mm2Ids(m2l)[components[[pid]]],cex = 1.0,pos=4)
				plot(m2l,id,title=titles[pid])
				layout(1)
			}else{
				plot(m2l,id,title=titles[pid])
			}
		}else{
			if(!is.null(components)){
				layout(matrix(c(1,2,2,2),byrow=TRUE,ncol=4))
				opar <- par("mar")
				par(mar=c(0.1,0.1,2.1,0.1))
				plot(0,xaxt="n",yaxt="n",xlab="",ylab="",ylim=c(0,1),xlim=c(0,1),type="n",bty="n")
				text(x = rep(0,length(components[[pid]])),y=seq(0.1,0.9,length=length(components[[pid]])),
					 labels = mm2Ids(m2l)[components[[pid]]],cex = 1.0,pos=4)
				plot(m2l,id)
				layout(1)
				par(mar=opar)
			}else{
				plot(m2l,id)
			}
		}

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
