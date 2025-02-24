#' @include MzDiffFormula.R


###CORRESPONDANCE TABLE
CORR_TABLE <- list("D"="dags","S"="spectra","P"="patterns","L"="losses")#"F"="fragments"
UNKNOWN_FORMULA <- NA_character_

###Accessors

#' Methods on ms2lib objects
#' 
#' Different atttributes from a ms2Lib object can be accessed using methods.
#' 
#' @rdname ms2Lib-methods
#' 
#' @export
#' @examples 
#' 
#' spectra <- mm2Spectra(m2l)
setMethod("mm2Spectra","ms2Lib",function(m2l){
	return(m2l@spectra)
})

#' @rdname ms2Lib-methods
#' 
#' @export
#' @examples 
#' 
#' infos <- mm2SpectraInfos(m2l)
setMethod("mm2SpectraInfos","ms2Lib",function(m2l){
	return(m2l@spectraInfo)
})

#' @rdname ms2Lib-methods
#' 
#' @export
#' @examples 
#' 
#' dags <- mm2Dags(m2l)
setMethod("mm2Dags","ms2Lib",function(m2l){
	return(m2l@dags)
})

#'@export
setMethod("mm2Ids","ms2Lib",function(m2l){
	return(m2l@ids)
})

#' @rdname ms2Lib-methods
#' 
#' @export
#' @examples 
#' 
#' atoms <- mm2Atoms(m2l)
setMethod("mm2Atoms","ms2Lib",function(m2l){
  return(m2l@atoms)
})

#' @rdname ms2Lib-methods
#' 
#' @export
#' @examples 
#' 
#' edge_labels <- mm2EdgesLabels(m2l)
setMethod("mm2EdgesLabels","ms2Lib",function(m2l){
	return(m2l@losses)
})

setMethod("mm2NodesLabels","ms2Lib",function(m2l){
	return(m2l@fragments)
})

#' @rdname ms2Lib-methods
#' 
#' @export
#' @examples 
#' 
#' patterns <- mm2Patterns(m2l)
setMethod("mm2Patterns","ms2Lib",function(m2l){
	return(m2l@patterns)
})

#' @rdname ms2Lib-methods
#' 
#' @export
#' @examples 
#' 
#' reduced_patterns <- mm2ReducedPatterns(m2l)
setMethod("mm2ReducedPatterns","ms2Lib",function(m2l){
	return(m2l@reducedPatterns)
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

setMethod("mm2Ids<-","ms2Lib",function(m2l, check=TRUE, value){
	if(check & any(startsWith(value,names(CORR_TABLE)))){
		stop("Forbidden prefixes for ids: ",paste(names(CORR_TABLE),collapse = ", "),
			 " found in ",value[startsWith(value,names(CORR_TABLE))])
	}

	if(length(value) != length(mm2Spectra(m2l))){
		stop("Number of furnished ids (",paste(length(value)),
			 ") should equal to spectra number (",paste(length(mm2Spectra(m2l))),")")
	}
	m2l@ids <- value
	###Check of the correctness of the IDs.
	m2l
})



#' @details For the `setIds` method, changing the ids of spectra can be useful to plot them.
#'
#' @param m2l An m2Lib object.
#' @param ids The new ids of an ms2Lib objects
#'
#' @export
#' 
#' @rdname ms2Lib-methods
#'
#' @examples
#' #Loading the data
#' data(m2l)
#' 
#' #Plot with all ids
#' plotOccurrences(m2l,"P15")
#' 
#' m2l <- setIds(m2l,paste("spectrum n",1:length(mm2Spectra(m2l)),sep=""))
#' 
#' #Plot with the new ids
#' plotOccurrences(m2l,"P15")
setMethod("setIds","ms2Lib",function(m2l,ids){
	mm2Ids(m2l) <- ids
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

setMethod("mm2ReducedPatterns<-","ms2Lib",function(m2l,value){
	m2l@reducedPatterns<- value
	m2l
})

setMethod("mm2Atoms<-","ms2Lib",function(m2l,value){
  m2l@atoms<- value
  m2l
})


isLoss <- function(m2l){
	m2l@loss
}


#' Get the supported input format for mineMS2
#' For now, the only supported format is mgf.
#' @return A character vector listing the supported format
#' @export
#'
#' @examples
#' recognisedFormat()
recognisedFormat <- function(){
	return(c("mgf"))
}


###Parse an mgf and return a list of the Spectrum2 objects
parse_mgf_spectrum2 <- function(filename){
	msnexp <- MSnbase::readMgfData(filename)
	lspec <- vector(mode="list",length=length(msnexp))
	for(i in 1:length(msnexp)){

		lspec[[i]] <- msnexp[[paste("X",i,sep="")]]
	}
	
	metadata <- fData(msnexp)[paste("X",1:length(msnexp),sep=""),]
	rlist <- list(spec=lspec,supp=metadata)

	return(rlist)
}

###Check the format of a list of files.
checkFormat <- function(x){
	rf <- recognisedFormat()
	splitted <- strsplit(x,".",fixed=TRUE)
	exts <- sapply(splitted,FUN = tail,1)
	is_ok <- exts %in% rf
	if(any(!(is_ok))) stop(paste("Unrecognized format:",exts[!is_ok],"in",x[!is_ok]))
	return(exts)
}

parseMS2file_line <- function(x){
	do.call(paste("parse",x[2],"spectrum2",sep="_"),list(filename=x[1]))
}

make_initial_title <- function(spec_infos){
	mz_str <- sprintf("%0.3f",spec_infos[,"mz.precursor"])
	paste("precursor mz: ",mz_str," (S",1:nrow(spec_infos),")",sep="")
}


convert_formula <- function(form_vec){
  ###We remove the white spaces and the special character
  form_vec <- trimws(str_replace(form_vec,"\\?|\\-|\\+|\\-",""))
  
  pnf <- which(is.na(form_vec)|(nchar(form_vec)==0))
  form_vec[pnf] <- UNKNOWN_FORMULA
  
  
  form_vec
}

#' ms2Lib constructor
#'
#' Create a ms2Lib object from files in the correct format. Eventually add some supplementary data given in a tsv file. 
#' This supplementary information may notably include the 'name' and the 'formula' of the molecules.
#'
#' @export
#' @param x May be one of the following
#' \itemize{
#' \item A character vector giving the path to a directory full of readable format.
#' \item A list of Spectrum2 objects which will be integrated directly.
#' \item A single .mgf file regrouping multiple spectra.
#' }
#' @param suppInfos Data frame containing the metadata to be associated to the spectra.
#' It should be of the same size as the number of spectra. If there is a 'file' column, this
#' column is used to match the file names. A 'formula' fields may also be present, and is then used to set the
#' formula of the molecules.
#' @param ids A supplementary vector giving a set of ids to design the spectra. It may be any character vector which
#' does not start with \code{'P', 'L', 'S'} as they are used internally by mineMS2. Alternatively if a suppInfos table is furnished
#' and if it contains an id fields, it will be used. If no ids are furnished, an id will be generated for each spectra in the form \code{'S1', 'S2', ..., 'SN'}
#' where N is the number of furnished spectra.
#' @param intThreshold The intensity threshold used to filter out the peaks.
#' @param infosFromFiles Shall the other information present in the files be added to the supplementary data (default: FALSE)
#' @aliases ms2Lib ms2Lib-constructor
#' @examples
#' #We locate the example file
#' path_demo <- system.file("dataset",package="mineMS2")
#' path_mgf <- file.path(path_demo,"dda_msms_pnordicum.mgf")
#' 
#' #Simple import
#' m2l <- ms2Lib(path_mgf)
#' 
#' #Import including some file formula
#' supp_infos_path <- file.path(path_demo,"dda_msms_pnordicum_supp.tsv")
#' supp_infos <- read.table(supp_infos_path,header=TRUE,sep=";")
#' m2l <- ms2Lib(path_mgf,suppInfos = supp_infos)
ms2Lib <- function(x, suppInfos = NULL,ids = NULL, intThreshold = NULL, infosFromFiles = FALSE){

	m2l <- new("ms2Lib")

	origin <- "R"
	lfiles <- NULL

	suppMetadata <- NULL
	###The kind of the acquisition is assessed there.
	if(is.list(x)){
		if(all(sapply(x,class) == "Spectrum2")){
			mm2Spectra(m2l) <- x
		}else{
			message("Unrecognized input, use one of the constructor described in the ms2Lib doc.")
		}

	}else if(is.character(x)){
		origin <- "file"
		if(length(x)==1){
		###Single mgf file or dir
			if(dir.exists(x)){
				lfiles <- list.files(x,full.names = TRUE) ## full path
				exts <- checkFormat(lfiles)
				message("Reading ",length(exts)," files with format(s): ",unique(exts))
				
				tres <- apply(matrix(c(lfiles,exts),ncol=2,byrow = FALSE),1,parseMS2file_line)
				mm2Spectra(m2l) <- sapply(tres,"[",i="spec")
				suppMetadata <- sapply(tres,"[",i="supp")
				
			}else{ ###Case of a single spectra.
				message("Reading MGF file ", x, ".")
				tres <- parseMS2file_line(c(x, 'mgf'))
				mm2Spectra(m2l) <- tres$spec
				suppMetadata <- tres$supp
				
			}
		}else{
			###Case of multiple individual spectra
			exts <- checkFormat(x)
			message("Reading ",length(exts)," files with format(s): ",unique(exts))
			
			mm2Spectra(m2l) <- apply(matrix(c(x,exts),ncol=2,byrow = FALSE),1,parseMS2file_line)
			
		}
	}

	mm2Spectra(m2l) <- do.call("c",mm2Spectra(m2l))
	
	#if(length(mm2Spectra(m2l))>2000){
	#  stop("At the moment it is impossible ot process more than 2,000 spectra at the same time.")
	#}

	###data.frame is initialized. With the mass of the precursors
	temp_df <- data.frame("mz.precursor" = sapply(mm2Spectra(m2l),function(x){
												  precursorMz(x)}))

	if(!is.null(intThreshold)){
		message("Removing peaks with an intensity lower than ",intThreshold)
		for(is in seq_along(mm2Spectra(m2l))){
			m2l@spectra[[is]] <- removePeaks(m2l@spectra[[is]],t = intThreshold)
			m2l@spectra[[is]] <- clean(m2l@spectra[[is]],all=TRUE)
		}
	}

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
	
	##Adding the supplementary information if necessary while check for an id fields.
	if(!is.null(suppInfos)){
		if(nrow(suppInfos)!= length(m2l@spectra)){
			stop("The number of suppInfos rows (",nrow(suppInfos),
				 ") do not match the number of spectra (",
				 length(m2l@spectra),")furnished")
		}else{
			## add a N column if not
			if(!("N" %in% colnames(suppInfos)))
			{	
				df_N <- data.frame(N = seq(1, nrow(suppInfos)))
				suppInfos <- cbind(suppInfos, df_N)
			}
			# ## if name or compounds column
			# if("name" %in% colnames(suppInfos))
			# {
			# 	colnames(suppInfos)[which(colnames(suppInfos) == "name")] <- "name"
			# }
			# if("compound" %in% colnames(suppInfos))
			# {
			# 	colnames(suppInfos)[which(colnames(suppInfos) == "compound")] <- "name"
			# }
			# if("Compound" %in% colnames(suppInfos))
			# {
			# 	colnames(suppInfos)[which(colnames(suppInfos) == "Compound")] <- "name"
			# }

			if("file" %in% colnames(suppInfos)){
				pm <- match(temp_df$file,suppInfos$file)
				if(any(is.na(pm))) stop('"file" column provided, but there was an error matching it against files.')

				temp_df <- cbind(temp_df,suppInfos[pm,])
			}else{
				temp_df <- cbind(temp_df,suppInfos)
			}
			if(("id" %in% colnames(suppInfos)) &
			   (is.null(ids))){
				mm2Ids(m2l) <- as.character(suppInfos[,"id"])
			}
		}
	}else{
		if(!is.null(ids)){
			mm2Ids(m2l) <- ids
		}else{
			mm2Ids(m2l,check=FALSE) <- paste("S",1:length(mm2Spectra(m2l)),sep="")
		}
	}
	
	if(infosFromFiles&!is.null(suppMetadata)){
	  temp_df <- cbind(temp_df,suppMetadata)
	}
	
	
	####Adding the molecular formula.
	cnames <- tolower(colnames(temp_df))
	
	pf <- which(cnames %in% c("formula","composition"))
	if(length(pf)==0){
	  message("No 'formula' column found. All formula are considered as unknown.")
	  temp_df$formula <- rep(UNKNOWN_FORMULA,nrow(temp_df))
	}else{
	  temp_df[,pf] <- convert_formula(as.character(temp_df[,pf]))
	}

	mm2SpectraInfos(m2l) <- temp_df
	m2l
}

#' Get all the precursor formulas of the spectra in the dataset
#' 
#' If given as supplementary information, this function returns all
#' the precursor formulas of the spectra in the dataset
#' 
#' @param m2l a ms2Lib object
#' 
#' @export
#' 
#' @examples 
#' data(m2l)
#' 
#' get_formula(m2l)
get_formula <- function(m2l){
  vf <- match("formula",tolower(trimws(colnames(m2l@spectraInfo))))
  return(m2l@spectraInfo[,vf])
}

getFormula <- function(m2l){
  vnames <- tolower(colnames(m2l@spectraInfo))
  return(as.character(m2l@spectraInfo[,match("formula",vnames)]))
}

#' Show an ms2Lib object.
#'
#' Get the string representation of an ms2Lib object
#'
#' @param object An m2Lib object to be shown.
#'
#' @return None.
#' @export
#'
#' @examples
#' #Loading the data
#' data(m2l)
#' m2l
setMethod("show","ms2Lib",function(object){
	cat("An ms2Lib object containing",length(object),"spectra.\n")
  if(length(mm2Atoms(object))!=0){
	  cat("It has",nrow(mm2EdgesLabels(object)),"edge labels built with atoms",paste(names(mm2Atoms(object)),collapse=","),".\n")
  }
  cat("The available spectral metadata are:",colnames(mm2SpectraInfos(object)),"\n")
	cat("It contains:",length(mm2Patterns(object)),"patterns\n")
	if(length(mm2ReducedPatterns(object))!=0) cat("It has been reduced to",
											   length(mm2ReducedPatterns(object)), "patterns")
})

#' Number of patterns of an ms2Lib object
#'
#' Return the number of mined patterns of an ms2Lib object
#'
#' @param x An m2Lib object to be shown.
#'
#' @return The number of mined patterns.
#' @export
#'
#' @examples
#' #Loading the data
#' data(m2l)
#' length(m2l)
setMethod("length","ms2Lib",function(x){
	return(length(mm2Spectra(x)))
})



#' Mine recurrent subgraphs from a set of graphs.
#'
#' Mine all the complete closed subgraphs from a set of preconstructed fragmentation graphs objects.
#'
#' @param m2l An ms2Lib object to be processed.
#' @param count The minimum number of spectra in which the spectrum need to be sampled.
#' @param sizeMin The minimum size of the mined patterns.
#' @param precursor Should only pattern coming from the root be mined.
#'
#' @return The filled ms2Lib object.
#' @export
#' 
#' @rdname mineClosedSubgraphs
#'
#' @examples
#' #Loading the data
#' data(m2l)
#' 
#' #Mining the subgraphs
#' m2l <- mineClosedSubgraphs(m2l,count=2,sizeMin = 1)
setMethod("mineClosedSubgraphs","ms2Lib",function(m2l, count = 2, sizeMin = 2, precursor = FALSE){
	if(count<2){
		warning("'count' parameter too low; value set to 2.")
		count <- 2
	}


	###Get the data.frame corresponding to the sizes.
	processing <- sapply(mm2Dags(m2l),function(x){
		ecount(x)>1
	})

	if(nrow(mm2EdgesLabels(m2l))==0){
		stop("No labels constructed, use the discretizeMzDifferences method first.")
	}


	kTree <- NULL
	if(nrow(mm2EdgesLabels(m2l))<600){
		kTree <- 2
	}else{
		kTree <- 1
	}

	if(sizeMin==1&nrow(mm2EdgesLabels(m2l))>600){
		###Wide variety of mass losses.
		warning("sizeMin parameters set to ", sizeMin, ": risk of computational overload.")
	}
	
	###We select the non empty graph to mine the patterns.
	sel_g <- which(sapply(mm2Dags(m2l),ecount)!=0)
	if(length(sel_g)==0) stop("No non-empty dags found.")

	###Converting the dags into data.frame.
	df_edges <- sapply(mm2Dags(m2l),fromIgraphToDf_edges,simplify = FALSE)[sel_g]
	df_vertices <-sapply(mm2Dags(m2l),fromIgraphToDf_vertices,simplify = FALSE)[sel_g]
	
	###Mining the patterns.
	resRcpp <- mineClosedDags(df_vertices,df_edges,processing,count,kTree,sizeMin,precursor)

	###Construction via fragPatternc constructor.
	mm2Patterns(m2l) <- sapply(resRcpp$patterns,function(x,sel_idx){
	  temp <- canonicalForm(fragPattern(x))
	  temp@occurrences[,1] <- sel_idx[temp@occurrences[,1]]
	  return(temp)
	 },USE.NAMES = FALSE,sel_idx=sel_g)

	###Initializing the names of the patterns.
	for(i in 1:length(m2l@patterns)) mm2Name(m2l@patterns[[i]]) <- paste("P",i,sep="")
	m2l@reducedPatterns <- seq_along(m2l@patterns)

	message("Processing completed: ", length(mm2Patterns(m2l)), " patterns mined.")
	m2l
})


###Parse an id
parseId <- function(m2l,idx){

	prefix <- substring(idx,1,1)
	number <- as.integer(substring(idx,2))
	if(!(prefix %in% names(CORR_TABLE))){
		###Checking if it's in the ids fields.
		resm <- match(idx,mm2Ids(m2l))
		if(is.na(resm)){
			stop("Invalid prefix ",prefix," authorized prefix are ",
			 	paste(names(CORR_TABLE),collapse=", "))
		}else{
			return(list(type=CORR_TABLE[["S"]],num=resm))
		}
	}

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
	if(is.data.frame(slot(m2l,arglist[[1]]))){
		return(slot(m2l,arglist[[1]])[arglist[[2]],])
	}
	(slot(m2l,arglist[[1]]))[[arglist[[2]]]]
}

#' Indexing function
#' 
#' Return an object given its index: can be a spectrum, a pattern, a DAG or a m/z difference.
#' 
#' @param x An ms2Lib object.
#' @param i The index, a string starting by S if it is a spectrum, P if it is a pattern D if it is a DAG, L if it is a m/z difference label
#' @param j unused.
#' @param drop unused.
#' @param ... unused.
#'
#' @return A list containing the object.
#' @export
#'
#' @examples
#' #Loading the data
#' data(m2l)
#' 
#' #Extracting a spectrum
#' m2l["S10"]
#' 
#' #The associated DAG
#' m2l["D10"]
#' 
#' #Extraction of a pattern
#' m2l["P54"]
#' 
#' #Extraction of m/z difference information
#' m2l["L20"]
setMethod('[','ms2Lib',function(x,i,j=NULL,...,drop=TRUE){
	if(length(i)==1){
		 temp <- mm2get(x,parseId(x,i))
		 # if(drop){
		 	return(temp)
		 # }
	}else{
		res <- lapply(lapply(i,parseId,m2l=x),mm2get,m2l=x)
		return(res)
	}
})


#' Plot an element given an idx
#'
#' The method depends of the type of the index furnished. The following prefixes are supported:
#' \code{'P'} for a pattern, the plot method of fragPattern object is called (using the package igraph).
#' \code{'D'} for a DAG, the plot_dag method is called using the package igraph.
#' \code{'S'} for a spectrum, a spectrum is plotted calling the Plot method of Spectrum2 object.
#' Any other value will be removed.
#'
#' @param x An ms2Lib oject.
#' @param y The index, a string starting by S if it a spectrum, P if it's a pattern or D if it's a dag.
#' @param title The title of the plot. If null, a title is automatically furnished.
#' @param tkplot if TRUE, for fragmentation graphs ("D"), the plot is dynamic and displayed using tkplot, otherwise it is static.
#' @param ... supplementary arguments to be passed by the method.
#'
#' @return a fragPattern object or nothing
#' @export
#'
#' @examples
#' #' #Loading the data
#' data(m2l)
#' 
#' #Plotting a pattern
#' plot(m2l,"S10")
#' 
#' #The associated DAG
#' plot(m2l,"D10")
#' 
#' #Plotting a pattern
#' plot(m2l,"P53")
setMethod("plot", "ms2Lib",
		  function(x,
		  		 y,title=NULL,tkplot=FALSE,
		  		 ...) {
		  	if(length(y)>1){
		  		warning("A single if may be plotted on each call, plotting the first element only")
		  		y <- y[1]
		  	}
		  	rid <- parseId(x,y)
		  	if(rid[[1]]=="patterns"){
		  	  toccs <- x[y]@occurrences[,1]
		  	  if(is.null(title)) title <- y
				return(plot(x[y],title = title,dags=mm2Dags(x),edgeLabels=(mm2EdgesLabels(x)),
				     atoms=names(x@atoms),formula=get_formula(x)[toccs],tkplot=tkplot,...))
		  	}else if(rid[[1]]=="spectra"){
		  		MSnbase:::plot_Spectrum2(x[y],full=TRUE,...)
		  	}else if(rid[[1]]=="dags"){
		  	  if(is.null(title)) title="Fragmentation Graph"
		  		plot_dag(x[y],idx=y,edgeLabels=(mm2EdgesLabels(x)),atoms=x@atoms,title=title,tkplot=tkplot,...)
		  	}else if(rid[[1]]=="losses"){
		  		stop("Impossible to plot a loss")
		  	}else if(rid[[1]]=="fragments"){
		  		stop("Impossible to plot a fragment")
		  	}
})

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

#' Finding a spectrum or a m/z difference
#' 
#' Search for a spectrum (S) or a m/z difference (L) in an ms2Lib object given a tolerance in ppm or dmz.
#'
#' @param m2l ms2Lib object
#' @param mz A double giving the mass to be searched.
#' @param type The data to be search, "S" for spectra and "L" for m/z differences
#' @param ppm The tolerance in ppm
#' @param dmz The minimum tolerance in Da, if the ppm tolerance is lower in Da than this threshold, this threshold is selected.
#'
#' @return A character vector giving the IDs of the found m/z differences or elements.
#' @export
#'
#' @examples
#' #' #Loading the data
#' data(m2l)
#' 
#' #Finding a spectrum (or several spectra) with a precusor mass of 391.1398 with a tolerance of 0.01
#' findMz(m2l, 391.1398, "S", dmz=0.01)
#' 
#' #Finding an m/z difference of 147 with a tolerance of 0.01
#' findMz(m2l, 147, "L",dmz=0.1)
findMz <- function(m2l,mz,type=c("S","L"),ppm=15,dmz=0.01){
	if(!is(m2l,"ms2Lib")) stop("m2l should be an 'ms2Lib' object.")
	type <- match.arg(type)
	tol <- max(dmz,mz*ppm*1e-6)
	if(type=="S"){
		if(length(mz)>1){
			return(sapply(mz,findMz.S,m2l=m2l,tol=tol,simplify=FALSE))
		}
		return(findMz.S(m2l,mz,tol))
	}else if(type == "L"){
		if(length(mz)>1){
			return(sapply(mz,findMz.L,m2l=m2l,tol=tol,simplify=FALSE))
		}
		return(findMz.L(m2l,mz,tol))
	}
	else{
		stop("The type of data to search must be a spectrum (S) or an m/z difference (L)")
	}
}


getInfo.L <- function(num,m2l){
	titles <- colnames(mm2EdgesLabels(m2l))
	titles <- titles[!(titles %in% 	c("sig", "fused",
									  "adv_loss", "pen_loss", "carb_only", "nitrogen_only", "ch_only",
									  "full_labels", "labs"))]
	return(mm2EdgesLabels(m2l)[
		num,titles,drop=FALSE])

}

getInfo.S <- function(num,m2l){
	titles <- colnames(mm2SpectraInfos(m2l))
	titles <- titles[!(titles %in% 	c("title"))]
	return(mm2SpectraInfos(m2l)[
		num,c(titles),drop=FALSE])
}

#' Available information 
#' 
#' Return the available information on an element of an ms2Lib object.
#'
#' @param m2l An ms2Lib object
#' @param ids A vector of IDs should be spectra or a m/z differences.
#' @param ... Supplementary information to be passed to the getInfo function.
#'
#' @return A data.frame giving information about the queried elements.
#' @export
#'
#' @examples
#' #Loading the data
#' data(m2l)
#' 
#' #Infos on 2 spectra
#' getInfo(m2l,c("S10","S25"))
#' 
#' #Infos on all the spectra
#' getInfo(m2l, "S")
#' 
#' #Infos on multiple losses
#' getInfo(m2l,c("L12","L42"))
#' 
#' #Example combined with the findMz 
#' getInfo(m2l,findMz(m2l,391.1398,dmz=0.01))
getInfo <- function(m2l,ids,...){
	if(!is(m2l,"ms2Lib")) stop("m2l should be an 'ms2Lib' object.")
	authorizedValue <- c("losses","spectra")

	if((length(ids)==1) && (nchar(ids)==1)){
		if(ids=="S"){
			ids <- vrange(m2l,"S")
		}else if(ids=="L"){
			ids <- vrange(m2l,"L")
		}else{
			stop("Single letters are only allowed for S and L.")
		}
	}

	pids <- sapply(ids,parseId,m2l=m2l,simplify=FALSE)

	type <- sapply(pids,'[[',i="type")

	if(all(type %in% authorizedValue)){
		num <- sapply(pids,'[[',i=2)
		res <- sapply(pids,function(x,m2l){
			if(x[[1]] == "losses"){
				return(getInfo.L(x[[2]],m2l))
			}
			if(x[[1]] == "spectra"){
				return(getInfo.S(x[[2]],m2l))
			}
		},simplify=FALSE,m2l=m2l)
		if(length(unique(type))==1){
			return(do.call("rbind",res))
		}else{
			return(res)
		}

	}else{
		stop("Invalid type for getInfo: ",unique(type[!(type %in% authorizedValue)]))
	}
}


#' Calculate the range of iterations
#' 
#' Return the full range of iterations for different objects for an MS2lib object.
#'
#' @param m2l An ms2Lib object
#' @param type "S" for spectra,"L" for m/z differences or "P" for patterns
#' @param reduced Used only if "type" is set to "P", shall the filtered pattern set be returned.
#' @param as.number If as.number is selected integer are returned without the prefix.
#'
#' @return A character vector giving the existing ids in the ms2Lib object,
#' @export
#'
#' @examples
#' #Loading the data
#' data(m2l)
#' 
#' #Range of iterations for spectra
#' vrange(m2l,"S")
setMethod("vrange","ms2Lib",function(m2l,type=c("S","L","P"), reduced=TRUE, as.number=FALSE){
	type <- match.arg(type)
	if(type=="S"){
		if(length(m2l)==0) return(character(0))
	  if(as.number) return(1:length(m2l))
		return(paste("S",1:length(m2l),sep=""))
	}
	if(type=="P"){
		if(length(mm2Patterns(m2l))==0) return(character(0))
		if(reduced){
		  if(as.number) return(mm2ReducedPatterns(m2l))
			return(paste("P",mm2ReducedPatterns(m2l),sep=""))
		}else{
		  if(as.number) return(1:length(mm2Patterns(m2l)))
			return(paste("P",1:length(mm2Patterns(m2l)),sep=""))
		}
	}
	if(type=="L"){
		if(nrow(mm2EdgesLabels(m2l))==0) return(character(0))
		return(paste("L",1:nrow(mm2EdgesLabels(m2l)),sep=""))
	}
	if(type=="F"){
		if(nrow(mm2NodesLabels(m2l))==0) return(character(0))
	  if(as.number) return(1:nrow(mm2NodesLabels(m2l)))
		return(paste("L",1:nrow(mm2NodesLabels(m2l)),sep=""))
	}
})


hasPatterns <- function(m2l){
  return(length(m2l@patterns)!=0)
}

setMethod("hasCoverage","ms2Lib",function(x){
  if(hasPatterns(x)){
    return(COVERAGE_NAME %in% colnames(x@patterns[[1]]@occurrences))
  }else{
    stop(paste("Patterns need to be computed before obtaining coverage.",sep=""))
  }
  return(FALSE)
})


#' Calculate the coverage
#' 
#' @param x The ms2Lib to be computed.
#'
#' @export
#' 
#' @rdname calculateCoverage
#'
#' @examples
#' #Loading the data
#' data(m2l)
#' 
#' #Calculate the coverage for an ms2Lib object
#' m2l <- calculateCoverage(m2l)
setMethod("calculateCoverage","ms2Lib",function(x){
  loss_mz <- mm2EdgesLabels(x)$mz
  mgs <- mm2Dags(x)
  pb <- txtProgressBar(min = 0, max = length(mm2Patterns(x)), initial = 0, char = "=",
                 width = NA, "Covergae calculation", "cov_calc", style = 3, file = "")
  for(i in seq_along(mm2Patterns(x))){
    setTxtProgressBar(pb, i)
    x@patterns[[i]] <- calculateCoverage(x@patterns[[i]],loss_mz,mgs)
  }
  return(x)
})
