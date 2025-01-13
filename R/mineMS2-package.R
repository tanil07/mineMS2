#' Mining of MS-MS spectra by FSM
#'
#' Mining of MS-MS spectra by frequent subgraph mining method.
#'
#' Mining of MS-MS spectra by frequent subgraph mining method followed. T
#'
#' "_PACKAGE"
#' @importFrom graphics plot layout text lines par points mtext
#' @importFrom grDevices rainbow
#' @importFrom MSnbase readMgfData precursorMz mz intensity removePeaks
#' @importMethodsFrom MSnbase plot clean fData
#' @importFrom stringr str_detect str_split str_sub fixed str_replace str_count
#' @importFrom Matrix Matrix
#' @importFrom utils tail txtProgressBar setTxtProgressBar
#' @importFrom stats density median na.omit pnorm
#' @importFrom igraph make_empty_graph graph.empty set_vertex_attr add_edges adjacent_vertices V E E<- V<- get.edge.ids edge_attr
#' edge_attr<- ecount as_data_frame graph_from_data_frame vcount plot.igraph tkplot layout_with_sugiyama layout_nicely vertex_attr %--%
#' head_of tail_of delete_edges set_edge_attr induced_subgraph incident incident_edges largest.cliques delete_vertices components write_graph
#' 
#' @import Rcpp
#' @importFrom ggplot2 ggplot geom_segment aes xlim ylim ggtitle theme scale_color_manual labs element_text 
#' @importFrom ggrepel geom_label_repel
#' @importFrom patchwork inset_element
#' @importFrom webchem get_csid cs_img
#' @importFrom png readPNG
#' @importFrom Spec2Annot find_compo_from_mass
#' @import grid
#' @importFrom gridExtra tableGrob
#' @name mineMS2-package
#' @useDynLib mineMS2
NULL