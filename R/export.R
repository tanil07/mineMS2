#' Export Loss Table
#'
#' @param m2l An ms2Lib object
#'
#' @return A data.frame with the losses
#' @export
#'
#' @examples
#' print("examples to be put here")
lossesTable <- function(m2l){
	if(nrow(mm2EdgesLabels(m2l))==0){
		stop("No loss labels found, use the 'discretizeMassLosses' loss functions")
	}
	return(mm2EdgesLabels(m2l)[,c("lab","mz","mzmin","mzmax","count","formula"),drop=FALSE])

}
