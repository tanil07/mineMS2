#' Export Losses Table
#'
#' Export the losses table aith the full informations
#' @param m2l An ms2Lib object
#'
#' @return A data.frame with the losses
#' @export
#'
#' @examples
#' #Loading the data
#' data(m2l)
#' 
#' head(lossesTable(m2l))
lossesTable <- function(m2l){
	if(nrow(mm2EdgesLabels(m2l))==0){
		stop("No loss labels found, use the 'discretizeMassLosses' loss functions")
	}
	return(mm2EdgesLabels(m2l)[,c("lab","mz","mzmin","mzmax","count","formula"),drop=FALSE])

}

#' Export Losses Table with information about the patterns containing those losses
#'
#' 
#' @param m2l An ms2Lib object
#'
#' @return A data.frame with the losses
#' @export
#'
lossesTableComplete <- function(m2l){
	if(nrow(mm2EdgesLabels(m2l))==0){
		stop("No loss labels found, use the 'discretizeMassLosses' loss functions")
	}
    if(length(mm2Patterns(m2l)) == 0)
    {
        stop("No patterns")
    }

    all_pertes <- mm2EdgesLabels(m2l)[,c("lab","mz","mzmin","mzmax","count","formula"),drop=FALSE]
	groups_patterns <- rep("", length(all_pertes[,"lab"]))
    nb_patterns <- rep(0, length(all_pertes[,"lab"]))

    patterns <- mm2Patterns(m2l) ## all patterns

    i = 1
    for(p in patterns)
    {
        print(i)
        g <- mm2Graph(p) ## graph of a pattern
        num_pertes <- unique(edge_attr(g, name="lab"))
        groups_patterns[num_pertes] <- paste(groups_patterns[num_pertes], i, sep=",")
        nb_patterns[num_pertes] <- nb_patterns[num_pertes] + 1
        i <- i+1
    }

    all_pertes[, "nb patterns"] <- nb_patterns
    all_pertes[, "patterns"] <- groups_patterns
    return(all_pertes)
}
