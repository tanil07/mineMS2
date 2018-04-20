###Motif plotting function uniquely.


###CONSTANTS
HIGH_MASS <- "high m/z"
MULTIPLE_FORMULA <- "2+ formula"


###Function which map a pattern on an occurences.
matchMzs <- function(mz1,mz2,ppm=15,dmz=0.03){
	o1 <- order(mz1)
	o2 <- order(mz2)
	closeMatch(mz1[o1],mz2[o2],o1,o2,ppm,dmz)
}


###Mapp the peaks of a pattern on another values.
get_mapping <- function(mg,patg,loss_mass,root=0,tol=0.005,ppm=10){

	###Raw matching by default
	mode <- match.arg(mode)

	###Extracting the correct root idx
	rootidx <- V(mg)[root+1]

	all_mz <- vertex_attr(mg,"mz")

	###We get the mzvalue of all the value after the root
	rootmz <- all_mz[root+1]
	lower_mz <- all_mz[(root+1):length(all_mz)]

	###We calculate the mass difference between all the nodes values.
	diff_mz <- rootmz-lower_mz

	###We alwyas match the closest value

	###Computing the root exts labels mg
	m0_idx <- mg[[rootidx,,edges=TRUE]][[1]]
	m0_mz <- edge_attr(mg,"lab",m0_idx)

	###Computing the root exts labels of pattern
	p0_idx <- patg[[1,,edges=TRUE]][[1]]
	p0_labs <- edge_attr(patg,"lab",index=p0_idx)
	p0_mz <- loss_mass[p0_labs]

	###Now we just have to match the value
	res_match <- matchMzs(lower_mz,p0_mz,50,0.1)

	###Normally there is a match for every value.
	###However if there is a mistake we remove it.
	nna <- which(!is.na(res_match))
	return(matrix(c(1,as.numeric(head_of(pat,p0_idx[nna])),
					as.numeric(rootidx),as.numeric(head_of(mg,m0_idx[res_match[nna]]))),
					byrow = TRUE,nrow=2))
}





# plot_pattern <- function(fp,edgeLabels,
# 						 nodeLab = c("default","label"),
# 						 edgeLab = c("none","formula"),
# 						 mzdigits = 3,vertex_size=30,tkplot=FALSE,...){
# 	nodeLab <- match.arg(nodeLab)
# 	edgeLab <- match.arg(edgeLab)
#
# 	g <- mm2Graph(fp)
#
# 	###We determine the edge label
# 	elabs <- edge_attr(g,"lab")
# 	txtlabs <- rep("",length(elabs))
# 	if(edgeLab=="formula")	txtlabs <- edgeLabels[elabs,"labs"]
#
# 	###Node labels, nodes are always in the right order
# 	nodes0 <- g[[1,]]
# 	labs0 <- edge_attr(g,"lab",g[[1,,edges=TRUE]][[1]])
# 	p0 <- 1
# 	pnon0 <- as.numeric(nodes0)
# 	nlabs <- rep("M",length(nodes0)+1)
# 	if(nodeLab=="label"){
# 		nlabs[pnon0] <- paste(nlabs[pnon0],labs0)
# 	}else if(nodeLab=="default"){
# 		nlabs[pnon0] <- edgeLabels[full_labs[labs0]]
# 	}
#
#
# 	###Title
# 	title <- paste("Pattern: ",mm2Name(fp))
#
# 	if(tkplot){
# 		tkplot(g,layout=(layout_with_sugiyama(g)$layout),canvas.width = 600,
# 			   canvas.height = 600,vertex.label=nlabs,
# 			   vertex.size=vertex_size,edge.label = txtlabs,
# 			   vertex.color="orange",...)
#
# 	}else{
# 		plot.igraph(g,layout=(layout_with_sugiyama(g)$layout),
# 					vertex.label=nlabs,vertex.size=vertex_size,edge.label = txtlabs,
# 					vertex.color="orange",...)
# 	}
# }



setMethod("plot", "fragPattern",
		  function(x,
		  		 y = NULL,
		  		 edgeLabels = NULL,
		  		 nodeLab = c("default", "label"),
		  		 edgeLab = c("formula","none"),
		  		 mzdigits = 3,
		  		 vertex_size = 30,
		  		 tkplot = FALSE,
		  		 ...) {
		  	nodeLab <- match.arg(nodeLab)
		  	edgeLab <- match.arg(edgeLab)

		  	g <- mm2Graph(x)
		  	###We determine the edge label
		  	elabs <- edge_attr(g, "lab")
		  	txtlabs <- rep("", length(elabs))
		  	if (edgeLab == "formula")
		  		txtlabs <- edgeLabels[elabs, "labs"]

		  	###Node labels, nodes are always in the right order
		  	nodes0 <- g[[1,]][[1]]
		  	labs0 <- edge_attr(g, "lab", g[[1, , edges = TRUE]][[1]])
		  	p0 <- 1
		  	pnon0 <- as.numeric(nodes0)
		  	nlabs <- rep("M", length(nodes0) + 1)
		  	if (nodeLab == "label") {
		  		nlabs[pnon0] <- paste(nlabs[pnon0], labs0)
		  	} else if (nodeLab == "default") {
		  		nlabs[pnon0] <- paste(nlabs[pnon0],edgeLabels[labs0,"full_labels"])
		  	}

		  	###Title
		  	title <- paste("Pattern: ", mm2Name(x))

		  	if (tkplot) {
		  		tkplot(
		  			g,
		  			layout = (layout_with_sugiyama(g)$layout),
		  			canvas.width = 600,
		  			canvas.height = 600,
		  			vertex.label = nlabs,
		  			vertex.size = vertex_size,
		  			edge.label = txtlabs,
		  			vertex.color = "orange",
		  			...
		  		)

		  	} else{
		  		plot.igraph(
		  			g,
		  			layout = (layout_with_sugiyama(g)$layout),
		  			vertex.label = nlabs,
		  			vertex.size = vertex_size,
		  			edge.label = txtlabs,
		  			vertex.color = "orange",
		  			...
		  		)
		  	}

		  })


setMethod("plotOccurences","ms2Lib",function(m2l,patternidx,...){
	plot(1)
})
# setMethod("plot","fragPattern",function(x,y=NULL,mode=c("fixed","interactive"),
# 										 labsNodes=c(),
# 										 labEdges=c("both","mz","formula"),...){
#
#
#
#
#
# })


