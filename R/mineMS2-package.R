#' Mining of MS-MS spectra by FSM
#'
#' Mining of MS-MS spectra by frequent subgraph mining method.
#'
#' Mining of MS-MS spectra by frequent subgraph mining method followed. T
#'
#' @docType package
#' @importFrom graphics plot layout text lines
#' @importFrom MSnbase readMgfData precursorMz mz intensity removePeaks
#' @importMethodsFrom MSnbase plot clean fData
#' @importFrom stringr str_detect str_split str_sub fixed str_replace
#' @importFrom Matrix Matrix
#' @importFrom igraph make_empty_graph graph.empty set_vertex_attr add_edges adjacent_vertices V E get.edge.ids edge_attr
#' edge_attr<- ecount as_data_frame graph_from_data_frame vcount plot.igraph tkplot layout_with_sugiyama vertex_attr
#' head_of tail_of delete_edges set_edge_attr induced_subgraph incident incident_edges
#' @import Rcpp
#' @name mineMS2-package
#' @useDynLib mineMS2
NULL
