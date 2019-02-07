###Motif plotting function uniquely.


###CONSTANTS
HIGH_MASS <- "high m/z"
MULTIPLE_FORMULA <- "2+ formula"


layoutMatrix <- function(n)
{
	if (n == 1) {
		return(matrix(c(1)))
	}
	if (n == 2) {
		return(matrix(c(1, 2), nrow = (2)))
	}
	if (n == 3) {
		return(matrix(c(1, 2, 3), nrow = (3)))
	}
	if (n == 4) {
		return(matrix(c(1, 2, 3, 4), nrow = (2), byrow = TRUE))
	}
	if (n == 5) {
		return(matrix(c(1, 2, 3, 4, 5, 6), nrow = (2), byrow = TRUE))
	}
	if (n == 6) {
		return(matrix(c(1, 2, 3, 4, 5, 6), nrow = (2), byrow = TRUE))
	}
	if (n == 7) {
		return(matrix(c(1, 2, 3, 4, 5, 6, 7, 8, 9), nrow = (3),
					  byrow = TRUE))
	}
	if (n == 8) {
		return(matrix(c(1, 2, 3, 4, 5, 6, 7, 8), nrow = (2),
					  byrow = TRUE))
	}
	if (n == 9) {
		return(matrix(c(1, 2, 3, 4, 5, 6, 7, 8, 9), nrow = (3),
					  byrow = TRUE))
	}
	if (n == 10) {
		return(matrix(c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12),
					  nrow = (4), byrow = TRUE))
	}
	if (n == 11) {
		return(matrix(c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12),
					  nrow = (4), byrow = TRUE))
	}
	if (n == 12) {
		return(matrix(c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12),
					  nrow = (4), byrow = TRUE))
	}
}


###Function which map a pattern on an occurences.
matchMzs <- function(mz1,mz2,ppm=15,dmz=0.03){
	o1 <- order(mz1)
	o2 <- order(mz2)
	closeMatch(mz1[o1],mz2[o2],o1,o2,ppm,dmz)
}


###Mapp the peaks of a pattern on another values.
get_mapping <- function(mg,patg,loss_mass,root=0,tol=0.02,ppm=20){

	###Extracting the correct root idx
	rootidx <- V(mg)[root+1]

	all_mz <- vertex_attr(mg,"mz")

	###We get the mzvalue of all the value after the root
	rootmz <- all_mz[root+1]

	###Computing the root exts labels mg
	m0_idx <- mg[[rootidx,,edges=TRUE]][[1]]
	m0_labs <- edge_attr(mg,"lab",m0_idx)
	m0_mz <- loss_mass[m0_labs]

	###Computing the root exts labels of pattern
	p0_idx <- patg[[1,,edges=TRUE]][[1]]
	p0_labs <- edge_attr(patg,"lab",index=p0_idx)
	p0_mz <- loss_mass[p0_labs]

	###Now we just have to match the value
	res_match <- matchMzs(m0_mz,p0_mz,ppm,tol)

	###Normally there is a match for every value.
	###However if there is a mistake we remove it.
	nna <- which(!is.na(res_match))
	return(matrix(c(1,as.numeric(head_of(patg,p0_idx[res_match[nna]])),
					as.numeric(rootidx),as.numeric(head_of(mg,m0_idx[nna]))),
					byrow = TRUE,nrow=2))
}





#' plotting a fargPattern object.
#'
#' plot a fragPattern associated fragmentation
#' dag using igraph capabilities.
#'
#' @param x
#' @param y = NULL,
#' @param edgeLabels = NULL,
#' @param nodeLab = c("default", "label"),
#' @param edgeLab = c("formula","none"),
#' @param mzdigits = 3,
#' @param vertex_size = 30,
#' @param tkplot = FALSE,
#' @param ...
#'
#' @return nothing
#' @export
#'
#' @examples
#' print("Exmaples ot be put here.")
setMethod("plot", "fragPattern",
		  function(x,
		  		 y = NULL,
		  		 title=NULL,
		  		 edgeLabels = NULL,
		  		 nodeLab = c("default", "label"),
		  		 edgeLab = c("formula","none"),
		  		 mzdigits = 3,
		  		 vertex_size = 55,
		  		 subNodes = NULL,
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

		  	if(!is.null(subNodes)){
		  		compN <- 1:length(nlabs)
		  		compN <- compN[-subNodes]
		  		nlabs[compN] <- ""

		  		se <- head_of(g,E(g))
		  		te <- tail_of(g,E(g))
		  		w_edges <- which((!(se %in% subNodes))&
		  								(!(te %in% subNodes)))
				if(length(w_edges)>0){
		  			txtlabs[w_edge] <- ""
				}
		  	}

		  	###Title
		  	if(is.null(title)){
		  	  title <- paste("Pattern: ", mm2Name(x))
		  	}

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
		  			main="title",
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
		  			main=title,
		  			...
		  		)
		  	}

		  })
#' PLotting the occurences of a fragPattern object inside
#' an ms2Lib object.
#'
#' @param m2l An ms2lib object.
#' @param pidx A pattern id or a frag_pattern object.
#' @param titles A vector of titles to be used. A default tile include id
#' and precursor will bne used by default
#' @param byPage The maximum number of spectra to be plotted by page.
#' @param ... supplementary function to be passed to theplot function.
#'
#' @return
#' @export
#'
#' @examples
setMethod("plotOccurences", "ms2Lib", function(m2l,
											   pidx,
											   titles = NULL,
											   byPage = 6,
											   subOccs=NULL,
											   ...) {
	###Verifying that a correct id have been queried.
	fp <- NULL
	if (class(pidx) == "fragPattern") {
		fp <- pidx
	} else if (class(pidx) == "character") {
		if (length(pidx) > 1)
			stop("Only a single index may be plotted every times.")
		fp <- m2l[pidx]
		if (class(fp) != "fragPattern") {
			stop("Wrong id, only a fragPattern object may by used by plotOccurences.")
		}
	}

	if (is.null(titles)){
		if(all(startsWith(mm2Ids(m2l),"S"))){
			titles <- mm2SpectraInfos(m2l)$title
		}else{
			titles <- paste(mm2Ids(m2l)," (S",1:length(mm2Spectra(m2l)),")",sep="")
		}
	}

	###Extracting the occurences,dags and the graph.
	occs <- mm2Occurences(fp)
	mgs <- mm2Dags(m2l)
	g <- mm2Graph(fp)

	###Mass loss are used for matching.
	loss_mz <- mm2EdgesLabels(m2l)$mz


	col_vec <- rainbow(nrow(occs))

	if(is.null(subOccs)){
		subOccs <- 1:nrow(occs)
	}else{
		if(all(subOccs>=1)&all(subOccs<=nrow(occs))){
			occs <- occs[subOccs,,drop=FALSE]
		}
	}


	###THE IDX IS HANDLED BY C++ AND SHOULD BE CORRECT.
	occs_gid <- occs[, 1]
	occs_pos <- occs[, 2]

	layout(layoutMatrix(min(byPage, length(occs_gid))))

	res_plot <- vector(mode="list",length=nrow(occs))
	for (i in 1:length(occs_gid)) {
		gid <- occs_gid[i]
		pos <- occs_pos[i]
		map <- get_mapping(mgs[[gid]], g, loss_mz, root = pos)


		###Plotting of the spectra
		intv <- vertex_attr(mgs[[gid]], "rel_int")
		mzs <- vertex_attr(mgs[[gid]], "mz")
		ids <- V(mgs[[gid]])

		###Peaks are split between matched and non matched.
		matched_peaks_idx <- match(map[2,], ids)
		non_matched_peaks_idx <- seq(1, length(mzs))[-matched_peaks_idx]
		
		col_seq <- rep("#000000FF",length(mzs))
		col_seq[matched_peaks_idx] <- col_vec[subOccs[i]]
		
		res_plot[[i]] <- data.frame(int=intv,mz=mzs,col=col_seq)
		
		plot(
			mzs[non_matched_peaks_idx],
			intv[non_matched_peaks_idx],
			type = "h",
			col = "#000000FF",
			lwd = 3,
			xlim = c(0, max(mzs) * 1.05),
			ylim = c(0, max(intv * 1.05)),
			main = titles[gid],
			xlab = "m/z",
			ylab = "intensity",...
		)

		###We add the matched peaks.
		points(mzs[matched_peaks_idx],
			   intv[matched_peaks_idx],
			   type = "h",
			   col =col_vec[subOccs[i]],
			   lwd = 3)
	}
	layout(1)
	
	return(invisible(res_plot))
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


