###CORRESPONDANCE TABLE
CORR_TABLE <- list("D"="dags","S"="spectra","P"="patterns","L"="losses")#"F"="fragments"



###Accessors
setMethod("mm2Spectra","ms2Lib",function(m2l){
	return(m2l@spectra)
})

setMethod("mm2SpectraInfos","ms2Lib",function(m2l){
	return(m2l@spectraInfo)
})

setMethod("mm2Dags","ms2Lib",function(m2l){
	return(m2l@dags)
})


setMethod("mm2EdgesLabels","ms2Lib",function(m2l){
	return(m2l@losses)
})

setMethod("mm2NodesLabels","ms2Lib",function(m2l){
	return(m2l@fragments)
})

setMethod("mm2Patterns","ms2Lib",function(m2l){
	return(m2l@patterns)
})

setMethod("mm2Lattice","ms2Lib",function(m2l){
	return(m2l@lattice)
})

setMethod("mm2LatticeIdxItems","ms2Lib",function(m2l){
	return(m2l@latticeIdxItems)
})

setMethod("mm2ReducedPatterns","ms2Lib",function(m2l){
	return(m2l@reducedPatterns)
})

setMethod("mm2ReducedLattice","ms2Lib",function(m2l){
	return(m2l@reducedLattice)
})


###Setter
setMethod("mm2Spectra<-","ms2Lib",function(m2l,value){
	m2l@spectra <- value
	m2l
})

setMethod("mm2SpectraInfos<-","ms2Lib",function(m2l,value){
	m2l@spectraInfo <- value
	m2l
})

setMethod("mm2Dags<-","ms2Lib",function(m2l,value){
	m2l@dags <- value
	m2l
})

setMethod("mm2EdgesLabels<-","ms2Lib",function(m2l,value){
	m2l@losses <- value
	m2l
})

setMethod("mm2NodesLabels<-","ms2Lib",function(m2l,value){
	m2l@fragments <- value
	m2l
})

setMethod("mm2Patterns<-","ms2Lib",function(m2l,value){
	m2l@patterns<- value
	m2l
})

setMethod("mm2Lattice<-","ms2Lib",function(m2l,value){
	m2l@lattice<- value
	m2l
})

setMethod("mm2LatticeIdxItems<-","ms2Lib",function(m2l,value){
	m2l@latticeIdxItems<- value
	m2l
})

setMethod("mm2ReducedPatterns<-","ms2Lib",function(m2l,value){
	m2l@reducedPatterns<- value
	m2l
})

setMethod("mm2ReducedLattice<-","ms2Lib",function(m2l,value){
	m2l@reducedLattice<- value
	m2l
})

###Function to get a spectrum.


###A list of the recognised format
recognisedFormat <- function(){
	return(c("mgf"))
}


###Parse an mgf and return a list of the Spectrum2 object
parse_mgf_spectrum2 <- function(filename){
	msnexp <- MSnbase::readMgfData(filename)
	lspec <- vector(mode="list",length=length(msnexp))
	for(i in length(msnexp)){
		lspec <- msnexp[[i]]
	}
	return(lspec)
}

###Check the format fo a list of files.
checkFormat <- function(x){
	rf <- recognisedFormat()
	splitted <- strsplit(x,".",fixed=TRUE)

	exts <- sapply(splitted,tail,1)

	is_ok <- exts %in% rf

	if(any(!(is_ok))) stop(paste("Unrecognized format:",exts[!is_ok],"in",x[!is_ok]))

	return(exts)
}

parseMS2file_line <- function(x){
	do.call(paste("parse",x[2],"spectrum2",sep="_"),list(filename=x[1]))
}

make_initial_title <- function(spec_infos){
	mz_str <- sprintf("%0.3f",spec_infos[,"mz.precursor"])
	paste("Precursor : ",mz_str," (S",1:nrow(spec_infos),")",sep="")
}



#' ms2Lib constructor
#'
#' Create a ms2Lib object from files in suppoerted input formats.
#'
#' @export
#' @param x May be one of the following
#' \itemize{
#' \item A character vector giving the path to a drictory full of readable format.
#' \item A list of spectrum2 object which will be integrated directly.
#' \item A single .mgfspectrum regrouping multiple files.
#' }
#' @param suppInfos Supplementary information to be associated to the spectra.
#' It should be of the same size as the number of spectra. If there is a "file" column, this
#' column is used to match the file.
#' @examples
#' print("examples to be put here")
ms2Lib <- function(x, suppInfos = NULL){

	m2l <- new("ms2Lib")

	origin <- "R"
	lfiles <- NULL

	###The kind of the acquisition is assessed there.
	if(class(x)=="list"){
		if(all(sapply(x,class) == "Spectrum2")){
			mm2Spectra(m2l) <- x
		}else{
			message("Unrecognized input, use one of the constructor described in the ms2Lib doc.")
		}

	}else if(class(x)=="character"){
		origin <- "file"
		if(length(x)==1){
		###Single mgf file or dir
			if(dir.exists(x)){
				lfiles <- list.files(x,full.names = TRUE)
				exts <- checkFormat(lfiles)
				message("Reading ",length(exts)," files with format(s): ",unique(exts))
				mm2Spectra(m2l) <- apply(matrix(c(lfiles,exts),ncol=2,byrow = FALSE),1,parseMS2file_line)

			}else{ ###Case of a single spectra.
				ext <- checkFormat(x)
				message("Reading ",length(exts)," files with format(s): ",unique(exts))
				mm2Spectra(m2l) <- parseMS2file_line(c(x,ext))
			}
		}else{
			###Case of multiples singles spectra
			exts <- checkFormat(x)
			message("Reading ",length(exts)," files with format(s): ",unique(exts))
			mm2Spectra(m2l) <- apply(matrix(c(x,exts),ncol=2,byrow = FALSE),1,parseMS2file_line)
		}
	}

	mm2Spectra(m2l) <- do.call("c",mm2Spectra(m2l))

	###data.frame is initialized. With the mass of the precusors
	temp_df <- data.frame("mz.precursor" = sapply(mm2Spectra(m2l),function(x){
												  precursorMz(x)}))

	if(origin=="file"){
		if(is.null(lfiles)){
			if(length(x)!=nrow(temp_df)){
				temp_df$file <- rep(x,nrow(temp_df))
			}else{
				temp_df$file <- x
			}
		}else{
			temp_df$file <- lfiles
		}
	}

	temp_df$title <- make_initial_title(temp_df)

	##Adding the supplementary information if necessary.
	if(!is.null(suppInfos)){
		if(nrow(suppInfos)!= length(m2l@spectra)){
			stop("The number of suppInfos rows (",nrow(suppInfos),
				 ") do not match the number of spectra (",
				 length(m2l@spectra),")furnished")
		}else{
			if("file" %in% colnames(suppInfos)){
				pm <- match(temp_df$file,suppInfos$file)
				if(any(is.na(pm))) stop('"file" column furnished, but there was an error matching it against files.')

				temp_df <- cbind(temp_df,suppInfos[pm,])
			}else{
				temp_df <- cbind(temp_df,suppInfos)
			}
		}

	}
	mm2SpectraInfos(m2l) <- temp_df
	m2l
}

#' @export
setMethod("show","ms2Lib",function(object){
	cat("An ms2Lib object containing",length(object),"spectra.\n")
	cat("It has",nrow(mm2EdgesLabels(object)),"edges labels.\n")
	cat("The available supplementary informations are:",colnames(mm2SpectraInfos(m2l)),"\n")
	cat("It contains: ",length(mm2Patterns(object)),"patterns")
})

#' @export
setMethod("length","ms2Lib",function(x){
	return(length(mm2Spectra(x)))
})





#' Mine recurrent subgraph from a set of graphs.
#'
#' Mine all the complete recurring subgraphs.
#'
#' @param m2l An m2Lib object to be processed.
#' @param num The number of spectra in which the spectrum need to be sampled.
#' @param sizeMin The minimum size of the mined patterns.
#' @param precursor Should only the occurences coming form the root be conserved.
#'
#' @return The filled ms2Lib object.
#' @export
#'
#' @examples
#' print("examples to be put here")
setMethod("mineClosedSubgraphs","ms2Lib",function(m2l, num = 2, sizeMin = 2, precursor = FALSE){
	if(num<2){
		warning("num parameters set to ",num," it is therefore set to 2.")
		num <- 2
	}


	###Get the data.frame correspoding to the sizes.
	processing <- sapply(mm2Dags(m2l),function(x){
		ecount(x)>1
	})

	if(nrow(mm2EdgesLabels(m2l))==0){
		stop("No labels constructed, use the DiscretizeMallLosses function first.")
	}
	kTree <- 2
	if(nrow(mm2EdgesLabels(m2l))<600){
		kTree <- 3
	}else{
		kTree <- 2
	}



	if(sizeMin==1&nrow(mm2EdgesLabels(m2l))>600){
		###Wide variety of mass losses.
		warning("sizeMin parameters set to ",sizeMin," risk of computational overhead.")
	}

	###Converting the dags into data.frame.
	df_edges <- sapply(mm2Dags(m2l),fromIgraphToDf_edges,simplify = FALSE)
	df_vertices <-sapply(mm2Dags(m2l),fromIgraphToDf_vertices,simplify = FALSE)

	###Mining the patterns.
	resRcpp <- mineClosedDags(df_vertices,df_edges,processing,num,kTree,sizeMin,precursor)

	mm2LatticeIdxItems(m2l) <- resRcpp$items

	###Construction via fragPatternc constructor.
	mm2Patterns(m2l) <- sapply(resRcpp$patterns,fragPattern,USE.NAMES = FALSE)

	###Initializing the names of the patterns.
	for(i in 1:length(m2l@patterns)) mm2Name(m2l@patterns[[i]]) <- paste("P",i,sep="")


	###Buillding of the lattice
	mm2Lattice(m2l) <- graph_from_data_frame(resRcpp$edges,directed=TRUE,resRcpp$nodes)

	###We add a label filed which give all the values of this.


	message("Processing finished, ",length(mm2Patterns(m2l))," patterns mined.")
	m2l
})


###Parse an id.
parseId <- function(m2l,idx){
	prefix <- substring(idx,1,1)
	number <- as.integer(substring(idx,2))
	# browser()
	if(!(prefix %in% names(CORR_TABLE))) stop("Invalid prefix ",prefix," authorized prefix are ",
											  paste(names(CORR_TABLE),collapse=", "))
	###The case of the L prfix is handled directly.
	if(prefix=="L"){
		if(nrow(mm2EdgesLabels(m2l))<number) stop("Invalid id for mass_losses: ",number,".")
		return(list(type=CORR_TABLE[[prefix]],num=number))
	}

	if( (number<=length(slot(m2l,CORR_TABLE[[prefix]])))&
		(number>=1)){
		return(list(type=CORR_TABLE[[prefix]],num=number))
	}else{
		stop("Invalid id for ",CORR_TABLE[[prefix]]," : ",number,".")
	}
}

mm2get <- function(m2l,arglist){
	if(class(slot(m2l,arglist[[1]]))=="data.frame"){
		return(slot(m2l,arglist[[1]])[arglist[[2]],])
	}
	(slot(m2l,arglist[[1]]))[[arglist[[2]]]]
}

#' Return the dag correspodning to a spectra or the modif.
#'
#' Indexing function.
#'
#' @param x An ms2Lib oject.
#' @param i The index, a string starting by S if it a spectrum, P if it's a pattern D if it's a dag, L if it's a loss
#' F if it's a fragment.
#'
#' @return a fragPattern object of an igraph graph object.
#' @export
#'
#' @examples
setMethod('[','ms2Lib',function(x,i,j=NULL,...,drop=TRUE){
	mm2get(x,parseId(x,i))
})


#' Plot an element given an idx
#'
#' The method depends of the tyep of the ID furnished. The following prefixes are supported :
#' \begin{itemize}
#' \item \textbf{P} A pattern is called plot method of fragPattern object.
#' \item \textbf{S} A spectrum is plotted calling the Plot method of spectrum 2 object.
#' \end{itemize}
#' Any other value will be removed.
#'
#' @param x An ms2Lib oject.
#' @param y The index, a string starting by S if it a spectrum, P if it's a pattern or D if it's a dag.
#' @param ... supplementary arguments to be passed by the method.
#'
#' @return a fragPattern object of an igraph graph object.
#' @export
#'
#' @examples
setMethod("plot", "ms2Lib",
		  function(x,
		  		 y,
		  		 ...) {
		  	rid <- parseId(x,y)
		  	if(rid[[1]]=="patterns"){
				plot(x[y],edgeLabels=(mm2EdgesLabels(x)),...)
		  	}else if(rid[[1]]=="spectra"){
		  		plot(x[y],full=TRUE,...)
		  	}else if(rid[[1]]=="dags"){
		  		stop("DAGS plotting not implemented at the moment.")
		  	}
})

#
# setMethod('[[','MSMSacquisition',function(x,i,j,...,drop=TRUE){


####Search and info functions
findMz.S <- function(m2l,mz,tol){
	infos <- mm2SpectraInfos(m2l)
	matched <- which(abs(infos$mz.precursor-mz)<tol)
	if(length(matched)==0) return(character(0))
	paste("S",matched,sep="")
}

findMz.L <- function(m2l,mz,tol){
	infos <- mm2EdgesLabels(m2l)
	matched <- which(abs(infos$mz-mz)<tol)
	if(length(matched)==0) return(character(0))
	paste("L",matched,sep="")
}

#' Search in an ms2Lib object?
#'
#' Search an ms2Lib object givena tolerance in ppm or dmz.
#'
#' @param m2l ms2Lib object
#' @param mz A double giving the mass to be searched.
#' @param type The data to be search, "S" for spectra and "L" for losses
#' @param ppm The tolerance in ppm
#' @param dmz The minimum tolerance in Da, if the ppm tolerance is lower in Da than this threshold, this threshold is selected.
#'
#' @return A character vector giving the IDs of the found losses or elements.
#' @export
#'
#' @examples
#' print("examples to be put here")
setMethod("findMz","ms2Lib",function(m2l,mz,type=c("S","L"),ppm=15,dmz=0.01){
	type <- match.arg(type)
	tol <- max(dmz,mz*ppm*1e-6)
	if(type=="S"){
		if(length(mz)>1){
			return(sapply(mz,findMz.S,m2l=m2l,tol=tol,simplify=FALSE))
		}
		return(findMz.S(m2l,mz,tol))
	}else{
		if(length(mz)>1){
			return(sapply(mz,findMz.L,m2l=m2l,tol=tol,simplify=FALSE))
		}
		return(findMz.L(m2l,mz,tol))
	}
})

# setMethod("select","ms2Lib",function(m2l,...){})

getInfo.L <- function(num){
	titles <- colnames(mm2EdgesLabels(m2l))
	titles <- titles[!(titles %in% 	c("sig", "fused",
									  "adv_loss", "pen_loss", "carb_only", "nitrogen_only", "ch_only",
									  "full_labels", "labs"))]




	return(mm2EdgesLabels(m2l)[
		num,titles,drop=FALSE])

}

getInfo.S <- function(num){
	titles <- colnames(mm2SpectraInfos(m2l))
	titles <- titles[!(titles %in% 	c("title"))]
	return(mm2SpectraInfos(m2l)[
		num,c(titles),drop=FALSE])
}

#' Return the available informationon a lost of a spectra.
#'
#' @param m2l An ms2Lib object
#' @param ids A vector of IDs
#' @param ... Supplementary information ot be passed to the getInfo function.
#'
#' @return A data.frame giving informations about the queried elements.
#' @export
#'
#' @examples
#' print("Examples to be put here")
setMethod("getInfo","ms2Lib", function(m2l,ids){
	authorizedValue <- c("losses","spectra")

	pids <- sapply(ids,parseId,m2l=m2l,simplify=FALSE)

	type <- sapply(pids,'[[',i=1)

	if(all(type %in% authorizedValue)){
		num <- sapply(pids,'[[',i=2)
		res <- sapply(pids,function(x){
			if(x[[1]] == "losses"){
				return(getInfo.L(x[[2]]))
			}
			if(x[[1]] == "spectra"){
				return(getInfo.S(x[[2]]))
			}
		},simplify=FALSE)
		if(length(unique(type))==1){
			return(do.call("rbind",res))
		}else{
			return(res)
		}

	}else{
		stop("Invalid type for getInfo: ",unique(type[!(type %in% authorizedValue)]))
	}
})
