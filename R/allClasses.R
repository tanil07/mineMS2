#'A class storing a set of MS/MS spectra, their associated graphs,
#'and some supplementary information.
#'
#' @export
#' @slot spectraInfo A data.frame including the spectra information at minima the "mz.precursor" field
#' and informations added by the user eventually.
#' @slot spectra A list of spectra stored under the form of Spectrum2 object.
#' @slot ids A chracter vector containing the ids of spectra.
#' @slot dags A list storing the set of graphs object corresponding to the MS-MS spectra.
#' @slot losses The discretized edge labels.
#' @slot fragments The discretized nodes labels.
#' @slot patterns A list storing the fragPattern objects.
#' @slot atoms A list of the atoms used to build the formula as well as their maximum number.
#' @slot loss A boolean indicating if the object is built with losses or fragments.
#' @slot reducedPatterns The ids of the reduced patterns.
#' @slot k The maximum depth of the constructed k path tree.
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
setClass("MzDiffFormula",
         slot = list(
           formula = "matrix"
         ),
         prototype = list(
           formula = matrix(NA_integer_,nrow=0,ncol=0)
         )
)