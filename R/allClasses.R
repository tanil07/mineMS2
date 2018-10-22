####An object which store a the discrete labels of an ms spectra and a list of spectra as graph.

#'A class storing a set of MS-MS spectra, their associated graphs,
#'and some supplemnetary inofmration.
#'
#' @export
#' @slot spectraInfo A data.frame including the spectra information at minima the "mz.precursor" field
#' and informations added by the user eventually.
#' @slot spectra A list of spectra sotred under the form of Spectrum2 object.
#' @slot ids A chracter vector containing the ids of spectra.
#' @slot dags A list storing the set of graphs object corresponding to the MS-MS spectra.
#' @slot edgesLabels The discretized edges labels.
#' @slot nodesLabels The discretized nodes labels.
#' @slot patterns A list storing the fragPattern objects.
#' @slot lattice The graph reprsentation of the lattice of the patterns.
#' @slot latticeIdxItems The ids of the nodes which are items in MS2process.
#' @slot ReducedPatterns The ids of the reduced patterns.
#' @slot ReducedLattice The graph of the reduced lattice.
#' @slot k The count of reduced patterns.
#' @aliases ms2Lib
#' @exportClass ms2Lib
setClass(
	"ms2Lib",
	slot = list(
		spectraInfo = "data.frame",
		spectra = "list",
		ids = "character",
		dags = "list",
		losses = "data.frame",
		fragments = "data.frame",
		patterns = "list",
		loss = "logical",
		lattice = "ANY",
		latticeIdxItems = "integer",
		reducedPatterns = "integer",
		reducedLattice = "ANY",
		k = "integer"
	),
	prototype = list(
		spectraInfo = data.frame(),
		spectra = list(),
		ids = character(),
		dags = list(),
		losses = data.frame(),
		fragments = data.frame(),
		patterns = list(),
		loss = TRUE,
		lattice = make_empty_graph(),
		latticeIdxItems = integer(),
		reducedPatterns = integer(),
		reducedLattice = make_empty_graph(),
		k = as.integer(1)

	)
)



setClass(
	"fragPattern",
	slot = list(
		graph = "ANY",
		occurences = "matrix",
		root = "integer",
		name = "character"
	),
	prototype = list(
		graph = make_empty_graph(),
		occurences = matrix(0,nrow=0,ncol=2),
		root = as.integer(1),
		name = character()
	)
)
