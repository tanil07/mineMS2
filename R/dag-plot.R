

plot_dag <- function(g,idx=NULL,
		  		 title="Fragmentation graph ",
		  		 edgeLabels = NULL,
		  		 edgeLab = c("formula","none"),
		  		 mzdigits = 3,
		  		 vertex_size = 55,
				 subnodes=NULL,
				 layoutm = NULL,
		  		 tkplot = FALSE,
				 def_col = "orange",
		  		 ...) {
		  	edgeLab <- match.arg(edgeLab)

		  	###We determine the edge label
		  	elabs <- edge_attr(g, "lab")
		  	txtlabs <- rep("", length(elabs))
		  	if (edgeLab == "formula")
		  		txtlabs <- as.character(edgeLabels[elabs, "labs"])

		  	###Node labels, nodes are always in the right order
		  	nodeslabs <- sprintf(paste("%0.",mzdigits,"f",sep=""),vertex_attr(g,"mz"))

		  	###Title
		  	if(is.null(title)){
		  		title <- paste("Fragmentation graph ",sep="")
		  		if(!is.null(idx)){
		  			tilte <- paste(title,idx)
		  		}
		  	}

		  	v_colvec <- rep(def_col,length(nodeslabs))
		  	ve_colvec <- rep("black",length(nodeslabs))
		  	lv_colvec <- rep("black",length(nodeslabs))
		  	e_colvec  <- rep("black",length(txtlabs))

		  	if(!is.null(subnodes)){
		  		compN <- 1:length(nodeslabs)
		  		compN <- compN[-subnodes]
		  		# nodeslabs[compN] <- ""
		  		v_colvec[compN] <- "white"
		  		ve_colvec[compN] <- "lightgrey"
		  		lv_colvec[compN] <- "lightgrey"
		  		# cat("supsup")

		  		se <- head_of(g,E(g))
		  		te <- tail_of(g,E(g))
		  		w_edges <- which((!(se %in% subnodes))|
		  						 	(!(te %in% subnodes)))
		  		if(length(w_edges)>0){
		  			txtlabs[w_edges] <- ""
		  			e_colvec[w_edges] <- "lightgrey"
		  		}
		  	}


		  	if(is.null(layoutm)){
		  		layoutm <- layout_with_sugiyama(g)$layout
		  	}

		  	if (tkplot) {
		  		tkplot(
		  			g,
		  			layout = (layoutm),
		  			canvas.width = 600,
		  			canvas.height = 600,
		  			vertex.label = nodeslabs,
		  			vertex.size = vertex_size,
		  			edge.label = txtlabs,
		  			vertex.frame.color = ve_colvec,
		  			vertex.color = v_colvec,
		  			vertex.label.color = lv_colvec,
		  			edge.color= e_colvec,
		  			main="title",
		  			...
		  		)

		  	} else{
		  		plot.igraph(
		  			g,
		  			layout = (layoutm),
		  			vertex.label = nodeslabs,
		  			vertex.size = vertex_size,
		  			edge.label = txtlabs,
		  			vertex.frame.color = ve_colvec,
		  			vertex.color = v_colvec,
		  			vertex.label.color = lv_colvec,
		  			edge.color= e_colvec,
		  			main=title,
		  			...
		  		)
		  	}

		  }
