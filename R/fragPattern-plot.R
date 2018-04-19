###Filed used for the plotting of motifs uniquely.



###Function which map a pattern on an occurences.

library(MS2discrete)

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



make_label_loss <- function(m2l){

	lab_edges <- mm2EdgesLabels(m2l)[,c("mz","formula")]

	formula <- info_tab$formula
	# MS2process:::RDBE(MS2process:::stringToFormula(formula[[5]]))

	str_formula <- str_split(formula,pattern = fixed("|"))

	str_formula <- sapply(str_formula,function(x){

		###If no label have been found.
		if(x=="NA") return("No formula")
		if(length(x)==1) return(x)
		if(length(x)>1) return(NA_character_)
	})

	###
	labs <- ifelse(is.na(str_formula),sprintf("%0.3f",info_tab$mz),str_formula)
	return(data.frame(mz=info_tab[,"mz"],labs=labs))
}


###TODO pass it in a class
plot_pattern <- function(pattern,loss_lab,
						 node_lab = c("diffmz","formula","both"),
						 edge_lab = c("none","formula"),
						 mzdigits = 3,vertex_size=30,...){
	node_lab <- match.arg(node_lab)
	edge_lab <- match.arg(edge_lab)


	if(class(pattern)!="igraph"){
		pattern <- read_graph(pattern,format = "graphml")
	}

	###We determine the edge label
	elabs <- edge_attr(pattern,"lab")
	txtlabs <- rep("",length(elabs))
	if(edge_lab=="formula") txtlabs <- loss_lab[elabs,"labs"]

	###Node labels.
	nlabs <- vertex_attr(pattern,"dist_prec")
	p0 <- which(nlabs==0)
	pnon0 <- which(nlabs!=0)
	if(node_lab=="diffmz"){
		nlabs[p0] <- "M"
		nlabs[pnon0] <- paste("M -\n",sprintf(paste("%0.",mzdigits,"f",sep=""),loss_lab[nlabs[pnon0],"mz"]))

	}else if(node_lab=="formula"){
		ostr <- rep("M",length(nlabs))
		diff_str <- rep("",length(nlabs))
		diff_str[pnon0] <- paste("-\n",loss_lab[nlabs[pnon0],"labs"],sep="")
	}else if(node_lab=="both"){

		ostr <- rep("M",length(nlabs))
		diff_str <- rep("",length(nlabs))
		diff_str[pnon0] <- paste("-\n",loss_lab[nlabs[pnon0],"labs"],sep="")
		tnlabs <- paste(ostr,diff_str)
		temp = rep("",length(nlabs))
		nlabs[pnon0] <- paste(tnlabs[pnon0],"\n(-",sprintf(paste("%0.",mzdigits,"f",sep=""),loss_lab[nlabs[pnon0],"mz"]),")",sep="")
		nlabs[p0] <- "M"
	}


	###Title
	title <- paste("Pattern: ",graph_attr(pattern,"norm"))



	plot.igraph(pattern,layout=(layout_with_sugiyama(pattern)$layout),
				vertex.label=nlabs,vertex.size=vertex_size,edge.label = txtlabs,
				vertex.color="orange",...)
}



#' Ploting the graph of a fragPattern object.
#'
#' Plot the graph of a fragPattern object.
#'
#' @param x a fragPattern object.
#' @param ... supplementary arguement to be passed ot plot.igraph.
#'
#' @return
#' @export
#'
#' @examples
setMethod("plot","fragPattern",function(x,y=NULL,edgeLabels=NULL,...){

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


