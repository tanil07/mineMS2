###Getter and setter
###ms2Lib object
setGeneric("mm2Spectra", function(m2l) standardGeneric("mm2Spectra"))

setGeneric("mm2Dags", function(m2l) standardGeneric("mm2Dags"))

setGeneric("mm2SpectraInfos", function(m2l) standardGeneric("mm2SpectraInfos"))

setGeneric("mm2EdgesLabels", function(m2l) standardGeneric("mm2EdgesLabels"))

setGeneric("mm2NodesLabels", function(m2l) standardGeneric("mm2NodesLabels"))

setGeneric("mm2Spectra<-", function(m2l,value) standardGeneric("mm2Spectra<-"))

setGeneric("mm2Dags<-", function(m2l,value) standardGeneric("mm2Dags<-"))

setGeneric("mm2SpectraInfos<-", function(m2l,value) standardGeneric("mm2SpectraInfos<-"))

setGeneric("mm2EdgesLabels<-", function(m2l,value) standardGeneric("mm2EdgesLabels<-"))

setGeneric("mm2NodesLabels<-", function(m2l,value) standardGeneric("mm2NodesLabels<-"))

###fragPattern object
setGeneric("mm2Graph", function(pat) standardGeneric("mm2Graph"))

setGeneric("mm2Occurences", function(pat) standardGeneric("mm2Occurences"))

setGeneric("mm2Root", function(pat) standardGeneric("mm2Root"))

setGeneric("mm2Graph<-", function(pat,value) standardGeneric("mm2Graph<-"))

setGeneric("mm2Occurences<-", function(pat,value) standardGeneric("mm2Occurences<-"))

setGeneric("mm2Root<-", function(pat,value) standardGeneric("mm2Root<-"))



###ms2Lib generic

setGeneric("discretizeMassLosses", function(m2l,...) standardGeneric("discretizeMassLosses"))

setGeneric("discretizeFragments", function(m2l,...) standardGeneric("discretizeFragments"))

###ms2Lib generic
setGeneric("mineClosedSubgraphs",function(m2l,...) standardGeneric("mineClosedSubgraphs"))
