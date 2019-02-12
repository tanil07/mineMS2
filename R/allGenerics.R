
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


###fragPattern object getter
setGeneric("mm2Graph", function(pat) standardGeneric("mm2Graph"))

setGeneric("mm2Occurences", function(pat) standardGeneric("mm2Occurences"))

setGeneric("mm2Root", function(pat) standardGeneric("mm2Root"))

setGeneric("mm2Name", function(pat) standardGeneric("mm2Name"))


###fragPattern object setter
setGeneric("mm2Graph<-", function(pat,value) standardGeneric("mm2Graph<-"))

setGeneric("mm2Occurences<-", function(pat,value) standardGeneric("mm2Occurences<-"))

setGeneric("mm2Root<-", function(pat,value) standardGeneric("mm2Root<-"))

setGeneric("mm2Name<-", function(pat,value) standardGeneric("mm2Name<-"))


###ms2Lib specific generics.
setGeneric("discretizeMassLosses", function(m2l,...) standardGeneric("discretizeMassLosses"))

setGeneric("discretizeMassFragments", function(m2l,...) standardGeneric("discretizeMassFragments"))

setGeneric("mineClosedSubgraphs",function(m2l,...) standardGeneric("mineClosedSubgraphs"))

setGeneric("mineClosedSubgraphsByComponent",function(m2l,...) standardGeneric("mineClosedSubgraphsByComponent"))

setGeneric("plotOccurences",function(m2l,pidx,...) standardGeneric("plotOccurences"))

setGeneric("exportLattice",function(m2l,filename,...) standardGeneric("exportLattice"))

setGeneric("score",function(m2l,...) standardGeneric("score"))

setGeneric("scoreLosses",function(m2l,...) standardGeneric("scoreLosses"))

setGeneric("scorePatterns",function(m2l,...) standardGeneric("scorePatterns"))

setGeneric("setIds",function(m2l,...) standardGeneric("setIds"))

setGeneric("addFormula",function(x,lf,...) standardGeneric("addFormula"))

setGeneric("isSubformula",function(x,lf,...) standardGeneric("isSubformula"))

setGeneric("isSubformulaIdx",function(x,lf,...) standardGeneric("isSubformulaIdx"))


###ms2Lib research generics.
setGeneric("vrange",function(m2l,type,...) standardGeneric("vrange"))

###ms2LibSplit generics.
setGeneric("setComponents",function(m2l,components,...) standardGeneric("setComponents"))


