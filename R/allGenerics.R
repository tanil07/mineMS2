
###ms2Lib object getter
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

###ms2Lib object setter
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


###fragPattern object getter
setGeneric("mm2Graph", function(pat) standardGeneric("mm2Graph"))

setGeneric("mm2Occurences", function(pat) standardGeneric("mm2Occurences"))

setGeneric("mm2Root", function(pat) standardGeneric("mm2Root"))

setGeneric("mm2Name", function(pat) standardGeneric("mm2Name"))

setGeneric("mm2CanonicalForm", function(pat) standardGeneric("mm2CanonicalForm"))


###fragPattern object setter
setGeneric("mm2Graph<-", function(pat,value) standardGeneric("mm2Graph<-"))

setGeneric("mm2Occurences<-", function(pat,value) standardGeneric("mm2Occurences<-"))

setGeneric("mm2Root<-", function(pat,value) standardGeneric("mm2Root<-"))

setGeneric("mm2Name<-", function(pat,value) standardGeneric("mm2Name<-"))

setGeneric("mm2CanonicalForm<-", function(pat,value) standardGeneric("mm2CanonicalForm<-"))


###ms2Lib specific generics.

#' Discretize mass differences.
#' 
#' @description
#' `discretizeMassLosses` discretize the mass differences iun MS-MS spectra.
#' 
#' @param m2l A ms2Lib object
#' @param ... Supplmentary arguments passed to the discretization method
#'
#' @details
#' This is a generic function: methods can be defined for it directly
#' or via the [Summary()] group generic. For this to work properly,
#' the arguments `...` should be unnamed, and dispatch is on the
#' first argument.
setGeneric("discretizeMassLosses", function(m2l,...) standardGeneric("discretizeMassLosses"))

setGeneric("discretizeMassFragments", function(m2l,...) standardGeneric("discretizeMassFragments"))

#' Extract frequent closed subgraphs.
#' 
#' `mineClosedSubgraphs` Extract frequent closed subgraphs from an ms2Lib object.
#' 
#' @param m2l A ms2Lib object
#' @param ... Supplmentary arguments passed to the mining method
#'
#' @details
#' This is a generic function: methods can be defined for it directly
#' or via the [Summary()] group generic. For this to work properly,
#' the arguments `...` should be unnamed, and dispatch is on the
#' first argument.
setGeneric("mineClosedSubgraphs",function(m2l,...) standardGeneric("mineClosedSubgraphs"))

setGeneric("mineClosedSubgraphsByComponent",function(m2l,...) standardGeneric("mineClosedSubgraphsByComponent"))

#' Plot occurences of fragPattern in library object.
#' 
#' @description
#' `plotOccurences` PLot a pattern on his associated fragmentation spectrum.
#'
#' @param m2l A ms2Lib object
#' @param pidx A pattern id to be plotted
#' @param ... Supplmentary arguments passed to the plotting method
#'
#' @details
#' This is a generic function: methods can be defined for it directly
#' or via the [Summary()] group generic. For this to work properly,
#' the arguments `...` should be unnamed, and dispatch is on the
#' first argument.
setGeneric("plotOccurences",function(m2l,pidx,...) standardGeneric("plotOccurences"))

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
#' @param ... Supplemetary arguments passed to the setIds method
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



###ms2Lib research generics.

#' Set the componetns of an ms2 library object.
#' 
#' @description
#' `vrange` return the range of the differnet elemnts of an ms2Lib object.
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

###ms2LibSplit generics.
#' Set the componetns of an ms2 library object.
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

###LossFormula generic
setGeneric("idxFormula",function(x,lf,...) standardGeneric("idxFormula"))
