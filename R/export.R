#'Export a lattice
#'
#'Export a lattice correspoding to an ms2Lib object to a file.
#'
#' @param m2l A ms2Lib object.
#' @param file A file connection ot a string giving th ename of a file.
#' @param format The format of export in a file.
#'
#' @return Nothing
#' @export
#'
#' @examples
#' print("Examples to be put here")
setMethod("exportLattice","ms2Lib",function(m2l,filename,format="graphml",...){
	if(vcount(mm2Lattice(m2l))==0) stop("The lattice is empty. Please mine graphs before exporting.")
	write_graph(mm2Lattice(m2l),filename,format=format)
})

#' export loss table
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
