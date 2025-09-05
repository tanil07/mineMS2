###Convert components furnished to spectral ids.
#convertComponent <- function(comp,ids){
#	if(is.numeric(comp)){
#		trad <- ids[comp]
#		if(any(is.na(trad))) stop("Invalid component elements: ",comp[is.na(trad)])
#	}else{
#		trad <- match(comp,ids)
#		if(any(is.na(trad))) stop("Invalid component elements: ",comp[is.na(trad)])
#		return(paste("S",trad,sep=""))
#	}
#}



#' Find patterns maximizing a score for each component
#' 
#' Find patterns maximizing the score given in metrics for each component.
#' Threshold may be used to indicate that no match has been found.
#' Several best patterns can be returned (argument top)
#'
#' @param m2l An ms2Lib object
#' @param components A list containing the components to be explained.
#' @param metric The metric or list of metrics used to match a component and the pattern occurrence.
#' @param ref_label A set of reference label to be used if the extracted component corresponds to an id over the position of the MS-MS spetra in the file
#' @param threshold The threshold used to discriminate the dataset
#' @param reduced Shall only the filtered patterns be considered
#' @param top the number of best patterns to return (by default 1)
#'
#' @return A list containing one data.frame by component including the following informations:
#'   \code{id} The id of the matched patterns.
#'   \code{f1} The F1-score between the pattern and the component.
#'   \code{precision} The precision between the pattern and the component.
#'   \code{recall} The recall between the pattern and the component.
#'   \code{miss} The miss rate between the pattern and the component.
#' 	 \code{size} The number of vertices in the graph pattern
#' @export
#'
#' @examples
#' #Loading the data
#' data(m2l)
#' data(molnet_df)
#' 
#' #Finding explainig component
#' molnet_igraph <- igraph::graph_from_data_frame(molnet_df$edges,
#' 				directed = FALSE, vertices = molnet_df$vertices)
#' components <- findGNPSComponents(molnet_igraph)
#' print(length(components))
#' fp <- findPatternsExplainingComponents(m2l,components)
#' print(fp[[1]])
#' #Plotting the first explaining pattern
#' plot(m2l,fp[[1]][1,"id"])
findPatternsExplainingComponents <- function(m2l,components,metric=c("f1"),ref_label=NULL,threshold=NA_real_,
                                             reduced=TRUE, top = 1){


	###We first convert the ids of the components to correct format.
  if(!is.null(ref_label)){
    components <- sapply(components,function(x,ref){
      ref[x]
    },ref=ref_label,simplify = FALSE)
  }
	###We determine the best format for each found pattern.
	bpat <- sapply(components,function(x,m2l,type,reduced,top){
		find.patterns.class(m2l,x,type = type,full = FALSE,reduced=reduced, top=top)
	},m2l=m2l,type=metric,simplify = FALSE,reduced=reduced, top=top)

	###Postprocessing
	return(bpat)
}
