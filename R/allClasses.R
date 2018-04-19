####An object which store a the discrete labels of an ms spectra and a list of spectra as graph.

#'A class storing a set of MS-MS spectra, their associated graphs,
#'and some supplemnetary inofmration.
#'
#' @export
#' @slot spectraInfo A data.frame including the spectra information at minima the "mz.precursor" field
#' and informations added by the user eventually.
#' @slot spectra A list of spectra sotred under the form of Spectrum2 object.
#' @slot dags A list storing the set of graphs object corresponding to the MS-MS spectra.
#' @slot edges.labels The discretized edges labels.
#' @slot nodes.labels The discretized nodes labels.
#' @aliases ms2Lib
#' @exportClass ms2Lib
setClass(
	"ms2Lib",
	slot = list(
		spectraInfo = "data.frame",
		spectra = "list",
		dags = "list",
		edgesLabels = "data.frame",
		nodesLabels = "data.frame"
	),
	prototype = list(
		spectraInfo = data.frame(),
		spectra = list(),
		dags = list(),
		edges.labels = data.frame(),
		noddes.labels = data.frame()
	)

)

setClass(
	"latticePatterns",

)


setClass(
	"fragPattern",
	slot = list(
		graph = "ANY",
		occurences = "matrix",
		root = "integer"
	),
	prototype = list(
		graph = make_empty_graph(),
		occurences = matrix(0,nrow=0,ncol=2),
		root = as.integer(1)
	)
)
