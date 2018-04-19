
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
	return(m2l@edgesLabels)
})

setMethod("mm2NodesLabels","ms2Lib",function(m2l){
	return(m2l@nodesLabels)
})

setMethod("mm2Patterns","ms2Lib",function(m2l){
	return(m2l@patterns)
})

setMethod("mm2Lattice","ms2Lib",function(m2l){
	return(m2l@lattice)
})

setMethod("mm2LatticeIdxItems","ms2Lib",function(m2l){
	return(m2l@lattice)
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
	m2l@edgesLabels <- value
	m2l
})

setMethod("mm2NodesLabels<-","ms2Lib",function(m2l,value){
	m2l@nodesLabels <- value
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
		warning("num parameters set to",num,"it is therefore set to 2.")
		num <- 2
	}

	###Get the data.frame correspoding to the sizes.
	processing <- sapply(mm2Dags(m2l),function(x){
		ecount(x)>1
	})

	kTree <- 2
	if(nrow(mm2EdgesLabels(m2l))<600){
		kTree <- 3
	}else{
		kTree <- 2
	}

	if(sizeMin==1&nrow(mm2EdgesLabels(m2l))>600){
		###Wide variety of mass losses.
		warning("sizeMin parameters set to",sizeMin,"risk of computational overhead.")
	}

	###Converting the dags into data.frame.
	df_edges <- sapply(mm2Dags(m2l),fromIgraphToDf_edges,simplify = FALSE)
	df_vertices <-sapply(mm2Dags(m2l),fromIgraphToDf_vertices,simplify = FALSE)

	###Mining the patterns.
	resRcpp <- mineClosedDags(df_vertices,df_edges,processing,num,kTree,sizeMin,precursor)

	mm2LatticeIdxItems(m2l) <- resRcpp$items

	###Construction via fragPatternc constructor.
	mm2Patterns(m2l) <- sapply(resRcpp$patterns,fragPattern,USE.NAMES = FALSE)

	###Buillding of the lattice
	mm2Lattice(m2l) <- graph_from_data_frame(resRcpp$edges,directed=TRUE,resRcpp$nodes)


	message("Processing finished, ",length(mm2Patterns(m2l))," patterns mined.")
	m2l
})
