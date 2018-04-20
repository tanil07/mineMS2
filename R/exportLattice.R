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
	if(vcount(mm2Lattice(m2l))) stop("The lattice is empty. Please mine graphs before exporting.")
	write_graph(mm2Lattice(m2l),file,format=format)
})
