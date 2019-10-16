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
#' @slot losses The discretized edges labels.
#' @slot fragments The discretized nodes labels.
#' @slot patterns A list storing the fragPattern objects.
#' @slot ReducedPatterns The ids of the reduced patterns.
#' @slot k The count of reduced patterns.
#' @aliases ms2Lib-class
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
		atoms = "list",
		patterns = "list",
		loss = "logical",
		reducedPatterns = "integer",
		# reducedLattice = "ANY",
		k = "integer"
	),
	prototype = list(
		spectraInfo = data.frame(),
		spectra = list(),
		ids = character(),
		dags = list(),
		losses = data.frame(),
		fragments = data.frame(),
		atoms = list(),
		patterns = list(),
		loss = TRUE,
		reducedPatterns = integer(),
		k = as.integer(1)
	)
)

setClass("ms2LibSplit",
		 representation = representation(components = "list",
		 								idxpcomponents = "numeric"),
		 contains="ms2Lib")


setClass(
	"fragPattern",
	slot = list(
		graph = "ANY",
		occurences = "matrix",
		root = "integer",
		name = "character",
		canonicalForm = "character"
	),
	prototype = list(
		graph = make_empty_graph(),
		occurences = matrix(0,nrow=0,ncol=2),
		root = as.integer(1),
		name = character(),
		canonicalForm = character()
	)
)

# Formula is a named numeric.
setClass("LossFormula",
         slot = list(
           formula = "matrix"
         ),
         prototype = list(
           formula = matrix(NA_integer_,nrow=0,ncol=0)
         )
)