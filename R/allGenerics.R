
###ms2Lib object getter
setGeneric("mm2Spectra", function(m2l) standardGeneric("mm2Spectra"))

setGeneric("mm2Dags", function(m2l) standardGeneric("mm2Dags"))

setGeneric("mm2SpectraInfos", function(m2l) standardGeneric("mm2SpectraInfos"))

setGeneric("mm2EdgesLabels", function(m2l) standardGeneric("mm2EdgesLabels"))

setGeneric("mm2NodesLabels", function(m2l) standardGeneric("mm2NodesLabels"))

setGeneric("mm2Patterns", function(m2l) standardGeneric("mm2Patterns"))

setGeneric("mm2Lattice", function(m2l) standardGeneric("mm2Lattice"))

setGeneric("mm2LatticeIdxItems", function(m2l) standardGeneric("mm2LatticeIdxItems"))


###ms2Lib object setter
setGeneric("mm2Spectra<-", function(m2l,value) standardGeneric("mm2Spectra<-"))

setGeneric("mm2Dags<-", function(m2l,value) standardGeneric("mm2Dags<-"))

setGeneric("mm2SpectraInfos<-", function(m2l,value) standardGeneric("mm2SpectraInfos<-"))

setGeneric("mm2EdgesLabels<-", function(m2l,value) standardGeneric("mm2EdgesLabels<-"))

setGeneric("mm2NodesLabels<-", function(m2l,value) standardGeneric("mm2NodesLabels<-"))

setGeneric("mm2Patterns<-", function(m2l,value) standardGeneric("mm2Patterns<-"))

setGeneric("mm2Lattice<-", function(m2l,value) standardGeneric("mm2Lattice<-"))

setGeneric("mm2LatticeIdxItems<-", function(m2l,value) standardGeneric("mm2LatticeIdxItems<-"))



###fragPattern object getter
setGeneric("mm2Graph", function(pat) standardGeneric("mm2Graph"))

setGeneric("mm2Occurences", function(pat) standardGeneric("mm2Occurences"))

setGeneric("mm2Root", function(pat) standardGeneric("mm2Root"))

###fragPattern object setter
setGeneric("mm2Graph<-", function(pat,value) standardGeneric("mm2Graph<-"))

setGeneric("mm2Occurences<-", function(pat,value) standardGeneric("mm2Occurences<-"))

setGeneric("mm2Root<-", function(pat,value) standardGeneric("mm2Root<-"))



###ms2Lib specific generic.
setGeneric("discretizeMassLosses", function(m2l,...) standardGeneric("discretizeMassLosses"))

setGeneric("discretizeFragments", function(m2l,...) standardGeneric("discretizeFragments"))

setGeneric("mineClosedSubgraphs",function(m2l,...) standardGeneric("mineClosedSubgraphs"))

setGeneric("plotOccurences",function(m2l,fp,...) standardGeneric("plotOccurences"))
