# ms2Lib object getter ----

setGeneric("mm2Spectra", function(m2l) standardGeneric("mm2Spectra"))

setGeneric("mm2Dags", function(m2l) standardGeneric("mm2Dags"))

setGeneric("mm2Ids", function(m2l) standardGeneric("mm2Ids"))

setGeneric("mm2SpectraInfos", function(m2l) standardGeneric("mm2SpectraInfos"))

setGeneric("mm2EdgesLabels", function(m2l) standardGeneric("mm2EdgesLabels"))

setGeneric("mm2NodesLabels", function(m2l) standardGeneric("mm2NodesLabels"))

setGeneric("mm2Patterns", function(m2l) standardGeneric("mm2Patterns"))

setGeneric("mm2Lattice", function(m2l) standardGeneric("mm2Lattice"))

setGeneric("mm2LatticeIdxItems", function(m2l) standardGeneric("mm2LatticeIdxItems"))

setGeneric("mm2ReducedPatterns", function(m2l) standardGeneric("mm2ReducedPatterns"))

setGeneric("mm2ReducedLattice", function(m2l) standardGeneric("mm2ReducedLattice"))

setGeneric("mm2Atoms", function(m2l) standardGeneric("mm2Atoms"))

# ms2Lib object setter ----

setGeneric("mm2Spectra<-", function(m2l,value) standardGeneric("mm2Spectra<-"))

setGeneric("mm2Dags<-", function(m2l,value) standardGeneric("mm2Dags<-"))

setGeneric("mm2Ids<-", function(m2l,value,...) standardGeneric("mm2Ids<-"))

setGeneric("mm2SpectraInfos<-", function(m2l,value) standardGeneric("mm2SpectraInfos<-"))

setGeneric("mm2EdgesLabels<-", function(m2l,value) standardGeneric("mm2EdgesLabels<-"))

setGeneric("mm2NodesLabels<-", function(m2l,value) standardGeneric("mm2NodesLabels<-"))

setGeneric("mm2Patterns<-", function(m2l,value) standardGeneric("mm2Patterns<-"))

setGeneric("mm2Lattice<-", function(m2l,value) standardGeneric("mm2Lattice<-"))

setGeneric("mm2LatticeIdxItems<-", function(m2l,value) standardGeneric("mm2LatticeIdxItems<-"))

setGeneric("mm2ReducedPatterns<-", function(m2l,value) standardGeneric("mm2ReducedPatterns<-"))

setGeneric("mm2ReducedLattice<-", function(m2l,value) standardGeneric("mm2ReducedLattice<-"))

setGeneric("mm2Atoms<-", function(m2l,value) standardGeneric("mm2Atoms<-"))


# fragPattern object getter ----

setGeneric("mm2Graph", function(pat) standardGeneric("mm2Graph"))

setGeneric("mm2Occurrences", function(pat) standardGeneric("mm2Occurrences"))

setGeneric("mm2Root", function(pat) standardGeneric("mm2Root"))

setGeneric("mm2Name", function(pat) standardGeneric("mm2Name"))

setGeneric("mm2CanonicalForm", function(pat) standardGeneric("mm2CanonicalForm"))


# fragPattern object setter ----

setGeneric("mm2Graph<-", function(pat,value) standardGeneric("mm2Graph<-"))

setGeneric("mm2Occurrences<-", function(pat,value) standardGeneric("mm2Occurrences<-"))

setGeneric("mm2Root<-", function(pat,value) standardGeneric("mm2Root<-"))

setGeneric("mm2Name<-", function(pat,value) standardGeneric("mm2Name<-"))

setGeneric("mm2CanonicalForm<-", function(pat,value) standardGeneric("mm2CanonicalForm<-"))


# ms2Lib specific generics ----

## discretizeMzDifferences ----

#' Discretize m/z differences.
#' 
#' @description
#' `discretizeMzDifferences` discretize the m/z differences in MS/MS spectra.
#' 
#' @param m2l A ms2Lib object
#' @param ... Supplementary arguments passed to the discretization method
#'
#' @details
#' This is a generic function: methods can be defined for it directly
#' or via the [Summary()] group generic. For this to work properly,
#' the arguments `...` should be unnamed, and dispatch is on the
#' first argument.
setGeneric("discretizeMzDifferences", function(m2l,...) standardGeneric("discretizeMzDifferences"))

setGeneric("discretizeMzFragments", function(m2l,...) standardGeneric("discretizeMzFragments"))

## mineClosedSubgraphs ----

#' Extract frequent closed subgraphs.
#' 
#' `mineClosedSubgraphs` Extract frequent closed subgraphs from an ms2Lib object.
#' 
#' @param m2l A ms2Lib object
#' @param ... Supplementary arguments passed to the mining method
#'
#' @details
#' This is a generic function: methods can be defined for it directly
#' or via the [Summary()] group generic. For this to work properly,
#' the arguments `...` should be unnamed, and dispatch is on the
#' first argument.
setGeneric("mineClosedSubgraphs",function(m2l,...) standardGeneric("mineClosedSubgraphs"))

setGeneric("mineClosedSubgraphsByComponent",function(m2l,...) standardGeneric("mineClosedSubgraphsByComponent"))

setGeneric("exportLattice",function(m2l,filename,...) standardGeneric("exportLattice"))

setGeneric("score",function(m2l,...) standardGeneric("score"))

setGeneric("scoreLosses",function(m2l,...) standardGeneric("scoreLosses"))

setGeneric("scorePatterns",function(m2l,...) standardGeneric("scorePatterns"))

#' Set the ids of a MS2lib object
#' 
#' @description
#' `setIds` set the ids of an ms2Lib objects
#' 
#' @param m2l An ms2Lib object
#' @param ... Supplementary arguments passed to the setIds method
#'
#' @details
#' This is a generic function: methods can be defined for it directly
#' or via the [Summary()] group generic. For this to work properly,
#' the arguments `...` should be unnamed, and dispatch is on the
#' first argument.
setGeneric("setIds",function(m2l,...) standardGeneric("setIds"))

setGeneric("addFormula",function(x,lf,...) standardGeneric("addFormula"))

setGeneric("isSubformula",function(x,lf,...) standardGeneric("isSubformula"))

setGeneric("isSubformulaIdx",function(x,lf,...) standardGeneric("isSubformulaIdx"))

#' Calculate the relative intensity covered by patterns.
#' 
#' @description
#' `calculateCoverage` calculate the relative intensity covered by a pattern or a library of patterns.
#' 
#' @param x An ms2Lib of a fragPattern object
#' @param ... supplementary argument passed to the methods
#'
#' @details
#' This is a generic function: methods can be defined for it directly
#' or via the [Summary()] group generic. For this to work properly,
#' the arguments `...` should be unnamed, and dispatch is on the
#' first argument.
setGeneric("calculateCoverage",function(x,...) standardGeneric("calculateCoverage"))

setGeneric("hasCoverage",function(x,...) standardGeneric("hasCoverage"))

# ms2Lib research generics ----

#' Set the components of an ms2 library object.
#' 
#' @description
#' `vrange` return the range of the different elements of an ms2Lib object.
#'
#' @param m2l An ms2Lib object
#' @param type The type of range to be returned
#' @param ... supplementary argument passed to the methods
#'
#' @details
#' This is a generic function: methods can be defined for it directly
#' or via the [Summary()] group generic. For this to work properly,
#' the arguments `...` should be unnamed, and dispatch is on the
#' first argument.
setGeneric("vrange",function(m2l,type,...) standardGeneric("vrange"))

# ms2LibSplit generics ----

#' Set the components of an ms2 library object.
#' 
#' @description
#' `setComponents` Set the components of an ms2Lib object.
#'
#' @param m2l An ms2Lib object
#' @param components The components given as a list of integers
#' @param ... supplementary argument passed to the methods
#'
#' @details
#' This is a generic function: methods can be defined for it directly
#' or via the [Summary()] group generic. For this to work properly,
#' the arguments `...` should be unnamed, and dispatch is on the
#' first argument.
setGeneric("setComponents",function(m2l,components,...) standardGeneric("setComponents"))

# MzDiffFormula generic
setGeneric("idxFormula",function(x,lf,...) standardGeneric("idxFormula"))

# Plots ----

## plotOccurrences ----

#' Plot occurrences of fragPattern in library object.
#' 
#' @description
#' `plotOccurrences` Plot a pattern on its associated fragmentation spectrum.
#'
#' @param m2l A ms2Lib object
#' @param pidx A pattern id to be plotted
#' @param ... Supplementary arguments passed to the plotting method
#'
#' @details
#' This is a generic function: methods can be defined for it directly
#' or via the [Summary()] group generic. For this to work properly,
#' the arguments `...` should be unnamed, and dispatch is on the
#' first argument.
#' @rdname plotOccurrences
#' @export
setGeneric("plotOccurrences",function(m2l,pidx,...) standardGeneric("plotOccurrences"))

## plotPatterns ----

#' Plot patterns and spectra occurrences of an ms2Lib object
#'
#' @description
#' `plotPatterns` plot the pattern graphs and supports (spectra occurrences) of an ms2Lib object.
#'
#' @param m2l An ms2Lib object
#' @param ids The ids to be plotted or NULL to plot all ids
#' @param components An id giving the component of each spectrum to be plotted in first page
#' @param occurrences Shall the occurrences be plotted as well as the pattern?
#' @param full Shall the full patterns be plotted?
#' @param byPage Maximum number of occurrences by page
#' @param titles A vector giving the titles of the MS/MS spectra
#' @param path_inchi name of a tabular file containing the inchi keys of the molecules in a column named "name"; if this table is available, the 2D structures corresponding to the spectra will be retrieved from ChemSpider (webchem package) and displayed in the plot along the spectra (default NULL)
#' @param ggplot.l (default TRUE) Should the ggplot version be displayed?
#' @param infos_col (for ggplot version) columns names from the spectra information to print
#' @param save_dir name of the directory to store the pdf files (default 'none' in case of display only)
#' @param ... supplementary argument passed to the methods
#' @return Nothing
#' @examples
#' #Loading the data
#' data(m2l)
#' 
#' #plotting the patterns
#' plotPatterns(m2l,c("P30","P51"))
#' @rdname plotPatterns
#' @export
setGeneric("plotPatterns", function(m2l,
                                    ids = NULL,
                                    components = NULL,
                                    occurrences = TRUE,
                                    full = FALSE,
                                    byPage = 9,
                                    titles = NULL,
                                    v_ggplot = TRUE, 
                                    path_inchi = NULL,
                                    ggplot.l = TRUE,
                                    infos_col = NULL,
                                    save_dir = "none",
                                    ...) standardGeneric("plotPatterns"))