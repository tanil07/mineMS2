

setMethod("mm2Graph", "fragPattern",function(pat){
	return(pat@graph)
})

setMethod("mm2Occurences", "fragPattern",function(pat){
	return(pat@occurences)
})

setMethod("mm2Root", "fragPattern",function(pat){
	return(pat@root)
})

setMethod("mm2Graph<-", "fragPattern", function(pat,value){
	pat@graph <- value
	pat
})

setMethod("mm2Occurences<-", "fragPattern", function(pat,value){
	pat@occurences <- value
	pat
})
setMethod("mm2Root<-","fragPattern", function(pat,value){
	pat@root <- value
	pat
})


#' @export
setMethod("show","fragPattern",function(object){
	cat("A fragPattern object with",length(mm2Graph(object))-1,"losses occuring ",
		nrow(mm2Occurences(object)),"times.\n")
})




