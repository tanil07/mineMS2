####Eventually first removing isolated points.
findAllCliques <- function(net,minSize = 3,vname="cluster index"){
	g <- induced_subgraph(net_gnps,V(net_gnps))

	big_clique <- largest.cliques(g)[[1]]

	pos_list <- 1
	list_cliques <- vector(mode="list",length=16)

	while(length(big_clique)>=minSize){

		if(pos_list>length(list_cliques)){
			list_cliques <- c(list_cliques,vector(mode="list",length=length(list_cliques)))
		}
		list_cliques[[pos_list]] <- vertex_attr(g,name=vname,index=big_clique)
		pos_list <- pos_list+1
		###We remove the vertices
		g <- delete_vertices(g,big_clique)
		big_clique <- largest.cliques(g)[[1]]
	}

	return(list_cliques[1:max(pos_list-1,1)])
}


findConnectedComponents <- function(net,minSize = 2,vname="cluster index",...){
	comp <- components(net,mode="weak")
	lelem <- sapply(comp$csize,numeric,simplify=FALSE)

	vec_idx <- rep(1,length(lelem))


	###We create a component list.
	Vlist <- V(net)
	for(i in seq_along(comp$membership)){
		lelem[[comp$membership[i]]][vec_idx[comp$membership[i]]] <- i
		vec_idx[comp$membership[i]] <- vec_idx[comp$membership[i]]+1
	}

	###We remove all the components with a size inferior tp mineSize
	lelem <- lelem[sapply(lelem,length)>=minSize]



	lapply(lelem,function(x,vidx,vname,net){vertex_attr(net,vname,index=vidx[x])},
		   vidx=Vlist,vname=vname,net=net)
}

# rcomp <- findConnectedComponents(net_gnps)

#' Find the annotable component of a GNPS network.
#'
#' Extract the annotated compoents of a network, at the moment, the extracted compoents are limited ot the clique.
#'
#' @param net A gnps network
#' @param minSize The minimum size of the detected clique.
#' @param pariThreshold For all the similarities, for each dataset,
#'
#' @return A list of the componets to be checked.
#' @export
#'
#' @examples
#' print("Examples to be put here")
findGNPSComponents <- function(net,minSize = 3,pairThreshold=0.9,vname="cluster index",eattr = "cosine_score"){

	###Finding all cliques.
	cliques <- findAllCliques(net,minSize = minSize,vname = vname)

	###Finding all connected components.
	connected_components <- findConnectedComponents(net,vname=vname,minSize=2)

	###We check for all  the connected components that any of the connected_components is a clique.

	###All the components are compared sequentiral.
	seqSizeClique <- sapply(cliques,length)
	seqSizeConnected <- sapply(connected_components,length)

	maxSize <-max(seqSizeClique)



	to_rm <- numeric(0)

	###We compare all the clique of similar size
	for(s in maxSize:1){
		pCliques <- which(seqSizeClique==s)
		if(length(pCliques)==0) next
		pConnected <- which(seqSizeConnected==s)
		if(length(pConnected)==0) next


		###Esle we doe the comparison.
		for(p1 in pCliques){
			for(p2 in pConnected){
				if(setequal(cliques[[p1]],connected_components[[p2]]))
					to_rm <- c(to_rm,p2)
			}
		}
	}

	###We filter the connected components list.

	pComp <- 1:length(connected_components)
	if(length(to_rm)!=0)	pComp <- pComp[-to_rm]

	if(length(pComp)!=0) connected_components <- connected_components[pComp]


	####Similarities passing the pairThreshold threshold are extracted.

	###Remvoing all the vertex from cliques and 2 paris connected components
	Vlist <- V(net)
	Vattr <- vertex_attr(net,vname)

	pairs <- connected_components[sapply(connected_components,length)==2]

	pcliques <- unique(do.call("c",c(cliques,pairs)))

	in_cliques <- match(Vattr,pcliques)
	out_cliques <- which(is.na(in_cliques))

	###All the edges are considered
	E_edges <-  E(net)[ Vlist %--% Vlist[out_cliques] ]
	E_edges <- E_edges[which(
		as.numeric(edge_attr(net,eattr,E_edges))>pairThreshold)]

	###Creating the components from the edges.
	highsim <- sapply(E_edges,function(x,net){
		c(head_of(net,x),tail_of(net,x))
	},net=net,simplify=FALSE)

	highsim <- sapply(highsim,function(x,net,vname){
		vertex_attr(net,vname,index=x)
	},net=net,vname=vname,simplify=FALSE)

	###Return the components.
	return(c(cliques,connected_components,highsim))
}


mapNetworkMS2lib <- function(net,m2l,vname="cluster index"){

	vidx <- V(net_gnps)

	vdescriptor<- vertex_attr(net,vname)

	vdescriptor <- as.numeric(vdescriptor)


	list(m2l_to_net=vidx[match(seq_along(m2l@spectra),vdescriptor)],net_to_m2l=match(vdescriptor,seq_along(m2l@spectra)))
}


componentsToDf <- function(comp,pID=NULL,maxVal = NULL,score = NULL){
	if(is.null(maxVal)){
		maxVal <- max(sapply(comp,max))
	}
	if(!is.null(score)&(length(score)!=length(comp)))stop("score should be a numeric vector of the same size than comp")

	###For each elements we get the component
	tfreq <- table(do.call("c",comp))
	mcomp <- max(as.numeric(tfreq))
	num_pat <- ifelse(length(pID)==0,0,mcomp+1)
	matval <- matrix(0,nrow=maxVal,ncol=mcomp+2)

	vscore <- NULL
	if(!is.null(score)){
		vscore <- rep(-1,maxVal)
	}

	matcomp <- NULL
	if(length(pID)!=0){
		matcomp <- matrix(NA_character_,nrow=maxVal,ncol=num_pat)
		colnames(matcomp) <- c(paste("PattMaxFComp",1:mcomp,sep=""),"ConsensusPatt")
	}

	matval[,1] <- 1:maxVal

	###Setting the ID of
	colnames(matval) <- c("ID",paste("Comp",1:mcomp,sep=""),"Consensus")


	###Now we fill the components with the values.
	vidx <- rep(2,max(as.numeric(names(tfreq))))

	for(ic in seq_along(comp)){
		for(c in comp[[ic]]){
			matval[c,vidx[c]] <- ic
			if(length(pID)!=0)
				matcomp[c,vidx[c]-1] <- pID[ic]
			vidx[c] <- vidx[c] + 1
			if(!is.null(score)) vscore[c] <- max(vscore[c],score[ic],na.rm=TRUE)
		}
	}


	###We add the last consensus column
	rkey <- apply(matval[,2:(mcomp+1)],1,paste,collapse="|")
	rkey <- as.numeric(as.factor(rkey))
	matval[,"Consensus"] <- rkey

	###If necessary we also build the consensus patterns.
	if(length(pID)!=0){
		rkey <- apply(matcomp[,1:mcomp],1,paste,collapse="|")
		matcomp[,"ConsensusPatt"] <- rkey
	}

	#####Building the resulting data.frame.
	df <- as.data.frame(matval)

	if(length(pID)!=0){
		df <- cbind(df,as.data.frame(matcomp))
	}

	if(!is.null(score)){
		df$score <- vscore
	}
	df
}



##Given a list of components, return for each vertices the index of the ocmpoents to which the vertex belong.


makeIdxTable <- function(components,maxVertices=70){

	res <- vector(mode="list",length=maxVertices)
	for(ir in seq_along(res)){
		res[[ir]] <- numeric(4)
	}
	seqidx <- rep(0,maxVertices)

	for(icomp in seq_along(components)){
		for(ic in components[[icomp]]){
			if((seqidx[ic]+1)>length(res[[ic]])){
				res[[ic]] <- c(res[[ic]],rep(NA_real_,length(res[[ic]])))
			}
			seqidx[ic] <- seqidx[ic]+1
			res[[ic]][seqidx[ic]] <- icomp

		}
	}

	for(ir in seq_along(res)){
		if(seqidx[ir]>=1){
		res[[ir]] <- res[[ir]][1:seqidx[ir]]
		}else{
			res[[ir]] <- numeric()
		}
	}
	return(res)
}


#' Title
#'
#' @param components jkkj
#' @param net kjkj
#' @param patterns jkkj
#' @param copy kjkj
#' @keepattr kjdffkdj
#' @param sep_infos kjkj
#'
#' @returnjkjk
#' @export
#'
#' @examples
#' print("Examples to be put here")
annotateNetwork <- function(components,net,patterns,copy=TRUE,
							keepattr = c("Compound_Name","SpectrumID","DefaultGroups",
										 "cluster index","precursor mass"),sep_infos=","){

	if(!is.null(keepattr)){
		if(copy) net <- induced_subgraph(net, V(net))
		attrnames <- vertex_attr_names(net)
		to_rm <- attrnames[is.na(match(attrnames,keepattr))]

		for(trm in to_rm){
		 net <- delete_vertex_attr(net, trm)
		}
	}

	DEFAULT <- "white"
	COLS_SEQ <- c("blue","green","green","orange","red","yellow","purple","cyan")

	fullseq <- rep(COLS_SEQ,ceiling(length(components)/length(COLS_SEQ)))


	resComp <- makeIdxTable(components,vcount(net))
	###Generating the col str.
	colstr <- sapply(resComp,function(x,colseq,defcol){
		if(length(x)>0){
		return(paste(colseq[x],collapse = ","))
		}else{
			return(defcol)
		}
	},defcol=DEFAULT,colseq=COLS_SEQ)

	legstr <- paste('stripechart: colorlist="',colstr,'"',sep="")
	ids <- sapply(patterns,function(x,sep){
		paste(x[,"id"],collapse=sep)
	},sep=sep_infos)


	rlist <- list()
	for(ne in colnames(patterns[[1]])){
		rlist[[ne]] <- sapply(patterns,function(x,sep,cname){
			paste(x[,cname],collapse=sep)
		},sep=sep_infos,cname= ne)
	}

	###We extract the maximum F1-score from this dataset.
	maxF1 <- sapply(patterns,function(x,sep){
		max(x[,"f1"],collapse=sep)
	},sep=sep_infos)

	str_components <- sapply(components,function(x,sep){
		paste(x,collapse = sep)
	},sep=sep_infos)

	temp_leg <- sapply(resComp,function(x,attstr,sep){
		if(length(x)==0) return("")
		return(paste(attstr[x],collapse = sep))
	},attstr=str_components,sep=sep_infos)
	net <- set_vertex_attr(graph = net,name = "components",value = temp_leg)

	###The consensus F1 score.
	temp_leg <- sapply(resComp,function(x,attstr,sep){
		if(length(x)==0) return("")
		return(paste(attstr[x],collapse = sep))
	},attstr=legstr,sep=sep_infos)
	net <- set_vertex_attr(graph = net,name = "colorComponents",value = temp_leg)

	###We now add the supplementary informations

	for(ne in names(rlist)){
		temp_leg <- sapply(resComp,function(x,attstr,sep){
			if(length(x)==0) return("")
			return(paste(attstr[x],collapse = sep))
		},attstr=rlist[[ne]],sep=sep_infos)
		net <- set_vertex_attr(graph = net,name = ne,value = temp_leg)
	}

	###The consensus F1 score.
	temp_leg <- sapply(resComp,function(x,attstr,sep){
		if(length(x)==0) return("")
		return(paste(attstr[x],collapse = sep))
	},attstr=maxF1,sep=sep_infos)
	net <- set_vertex_attr(graph = net,name = "maxF1",value = temp_leg)
	return(net)
}



###score is a vector of the score of the same size than component.
annotateNetworkWithComponents <- function(net,comps,vname="cluster index",score=NULL,patterns=NULL){
	seqIdx <- vertex_attr(net,name = vname)
	vmax <- max(seqIdx)

	rdf <- NULL
	if(!is.null(patterns)){
		rdf <- componentsToDf(comps,maxVal = vmax,score=score,pID=patterns)
	}else{
		rdf <- componentsToDf(comps,maxVal = vmax,score=score)
	}

	###We get the correct order
	om <- match(seqIdx,rdf[,1])

	if(any(is.na(om))){
		###We remove the missing vertices
		pok <- which(!is.na(om))
		rdf <- rdf[pok,]
	}



	###Adding attributes to the graph.
	cnames <- colnames(rdf)

	colsComp <- 2:(ncol(rdf)-1)
	colsPat <- numeric()

	if(!is.null(patterns)){
		sscore <- ifelse(is.null(score),0,1)
		fac <- (ncol(rdf)-sscore-1)/2-1
		colsComp <- 2:(fac)
		colsPat <- (fac+2):(2*fac+1)
	}


	for(i in colsComp){
		net <- set_vertex_attr(graph = net,name = cnames[i],value = as.integer(rdf[om,i]))
	}

	for(i in colsPat){
		net <- set_vertex_attr(graph = net,name = cnames[i],value = rdf[om,i])
	}
	net <- set_vertex_attr(graph = net,name = "Consensus",value = as.character(rdf[om,"Consensus"]))
	if(!is.null(patterns)) net <- set_vertex_attr(graph = net,name = "ConsensusPatt",value = as.character(rdf[om,"ConsensusPatt"]))
	if(!is.null(score)) 	net <- set_vertex_attr(graph = net,name = "F1",value = as.numeric(rdf[om,"score"]))
	net
}






#' Find a pattern for each MS-MS network components.
#'
#' At the moment only clique are exported.
#'
#' @param net An MS-MS network
#' @param type The matching algorithm to be used to determine the best pattern
#' @param m2l An ms2Lib object
#' @param vnames The column used to determine the id.
#' @param ... Supplementary parameters to be passed to findGNPSComponents.
#'
#' @return
#' @export
#'
#' @examples
#' print("Example to be put here")
ExplainGNPSComponents<-function(net,m2l,metric=c("f1","recall","accuracy"),vname="cluster index",threshold=NA_real_){

	type <- match.arg(type)

	message("Extracting network components.")

	##Each compoent is retruned as a list of vname attributes value.
	rcomp <- findGNPSComponents(net,...)

	###For each components find the best matching pattern.
	message("Calculating best matching patterns.")


	bpat <- findPatternsFromComponents(m2l,rcomp,metric=metric,threshold=threshold)

	vscore <-  sapply(bpat,function(x){
		unique(x[,metric])
	})

	###We extract the pattern IDs if there is multiple ids.
	idpats <- sapply(bpat,function(x){
		paste(x[,"id"],sep="|",collapse = "|")
	})


	message("Annotating network.")
	net <- annotateNetworkWithComponents(net,rcomp,score=vscore,vname="cluster index",patterns = idpats)

	# bpat <- do.call("rbind",bpat)
	list(network = net,components_ms2lib = rpattern,components_graph=rcomp,score=vscore,
		 patterns = bpat)
}

