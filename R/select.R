####Implementation of the selection functions.


checkSize <- function(vec,pos){
  if(typeof(vec) %in% c("numeric","integer","character","logical")){
  	if(length(vec)<pos){
  		vec <- c(vec,do.call(typeof(vec),args=list(length(vec))))
  	}
  	return(vec)
  }
  
  if(is.data.frame(vec)){
    if(nrow(vec)<pos){
      vec <- vec[rep(1:(nrow(vec)),times=2),]
    }
    return(vec)
  }
  
  if(is.list(vec)){
    if(length(vec)<pos){
      vec <- c(vec,vector(mode="list",length=length(vec)))
    }
    return(vec)
  }
  
  stop("Unrecognized data type :",typeof(vec))
}

#' Selection function
#'
#' Return spectra objects, that contain particular patterns, or m/z differences objects that are present in particular patterns or spectra.
#'
#' @param m2l An ms2Lib object.
#' @param ids Valid IDs of objects (spectra or m/z differences)
#' @param vals the type of values to be returned. Only certain combinations are allowed: we can search for spectra (S) in patterns (P), or m/z differences (L) in spectra (S) or patterns (P)
#'
#' @return The objects of the corresponding type.
#' @export
#'
#' @examples
#' data(m2l)
#' 
#' ## Search for the patterns in which appear the spectrum S1
#' select(m2l, "S1", "P")
#' 
#' ## Search for the patterns in which appear the m/z difference L1
#' select(m2l, "L1", "P")
#' 
#' ## Search for the spectra in which appear the m/z difference L1
#' select(m2l, "L1", "S")
#' 
#' ## can be combined with the findMz function
#' select(m2l, findMz(m2l, mz=147, "L", ppm=8, dmz=0.3), "P")
select <- function(m2l,ids,vals=c("P","S")){
	if(!is(m2l,"ms2Lib")) stop("m2l should be an ms2Lib object.")

	vals <- match.arg(vals)

	ids <- sapply(ids,parseId,m2l=m2l,simplify = FALSE)

	type <- unique(sapply(ids,'[[',i=1))

	if(length(type)!=1) stop("All the ids should be of the same type.")

	nums <- sapply(ids,'[[',i=2)

	if(type=="spectra"){
		if(vals=="P"){
			return(sapply(nums,select.patterns.spectra,m2l=m2l,simplify = FALSE))
		}else{
			stop("not implemented")
		}
	}else if(type=="losses"){
		if(vals=="P"){
			return(sapply(nums,select.patterns.losses,m2l=m2l,simplify = FALSE))
		}else if(vals=="S"){
			return(sapply(nums,select.spectra.losses,m2l=m2l,simplify = FALSE))
		}else{
			stop("not implemented")
		}
	}else if(type=="patterns"){
		stop("not implemented")
	}else if(type=="fragments"){
		stop("not implemented")
	}
}

#' Find the pattern with the best coverage given a set of spectra.
#'
#' Find the pattern with the best coverage, which means that covered the most intensity, for a set of spectra.
#'
#' @param m2l an ms2Lib Object.
#' @param ids A list of ids to be searched for (for now, only spectra are implemented).
#' @param vals The type of the returned values (for now, only patterns are implemented (P)).
#' @param ... supplementary arguments.
#'
#' @return The ids of the patterns or spectra containing all the ids.
#' @export
#'
#' @examples
#' data(m2l)
#' ## The calculation of coverage must have been done before searching
#' m2l <- calculateCoverage(m2l)
#' find(m2l, "S7", "P")
find <- function(m2l,ids,vals=c("P","S"),...){
	if(!is(m2l,"ms2Lib")) stop("m2l should be an ms2Lib object.")

	vals <- match.arg(vals)

	ids <- sapply(ids,parseId,m2l=m2l,simplify = FALSE)

	type <- unique(sapply(ids,'[[',i=1))

	if(length(type)!=1) stop("All the ids should be of the same type.")
	if(!(type %in% c("patterns","spectra"))){
		stop("ids should correspond to spectra or patterns.")
	}
	nums <- unique(sapply(ids,'[[',i=2))
	if((type=="spectra") && (vals=="P")){
		return(find.patterns.spectra(m2l,nums,...))
	}

	if((type=="patterns") && (vals=="S")){
		stop("Not implemented (soon)")
	}
	stop("Comparison between ",type," and ",vals," not implemented.")
}


find.patterns.spectra<- function(m2l,ids,reduced=FALSE,metric = c("coverage","size","none"),partial = TRUE){
  
  metric <- match.arg(metric)
  
	idsp <- vrange(m2l,"P",reduced)
	
	if((!hasCoverage(m2l@patterns[[1]])) & (metric=="coverage")) stop("use calculateCoverage function before searching with type to coverage") 
	
	info_patt_full <- data.frame(id=character(10),size=numeric(10),inter=numeric(10),stringsAsFactors = FALSE)
	info_patt_partial <- data.frame(id=character(10),size=numeric(10),inter=numeric(10),stringsAsFactors = FALSE)
	
	if(hasCoverage(m2l@patterns[[1]])){
	  info_patt_full$coverage <- numeric(10)
	  info_patt_partial$coverage <- numeric(10)
	}

	nsol <- 0
	nsol_partial <- 0
	for(ip in seq_along(idsp)){
		temp_pat <- m2l[idsp[ip]]
		occs <- mm2Occurrences(temp_pat)
		vm <- (ids %in% occs[,"gid"])

		if(all(vm)){
			nsol <- nsol+1
			info_patt_full <- checkSize(info_patt_full,nsol)
			info_patt_full$id[nsol] <- idsp[ip]
			info_patt_full$size[nsol] <- vcount(temp_pat@graph)-1
			info_patt_full$inter[nsol] <- sum(vm)
			if(hasCoverage(m2l@patterns[[1]])){
			  info_patt_full$coverage[nsol] <- median(occs[,"coverage"])
			}
		}
		if(any(vm)){
		  nsol_partial <- nsol_partial+1
		  info_patt_partial <- checkSize(info_patt_partial,nsol_partial)
		  info_patt_partial$id[nsol_partial] <- idsp[ip]
		  info_patt_partial$size[nsol_partial] <- vcount(temp_pat@graph)-1
		  info_patt_partial$inter[nsol_partial] <- sum(vm)
		  
		  if(hasCoverage(m2l@patterns[[1]])){
		    info_patt_partial$coverage[nsol_partial] <- median(occs[,"coverage"])
		  }
		}
	}
	resids <- NULL
	if(nsol!=0){
	  info_patt_full <- info_patt_full[1:nsol,]
		if(metric != "none"){
		  resids <- info_patt_full$id[which.max(info_patt_full[[metric]])]
		}
	  if(metric == "none") resids <- info_patt_full$id[info_patt_full$inter>=2]
	}else{
		if(partial){
			if(nsol_partial != 0){
				resids <- info_patt_partial$id[which.max(info_patt_full[[metric]])]
			}
		}
	}
	
	return(resids)
}


#' List all patterns ("P") found for each spectrum (S) or loss (l)
#'
#' List all the patterns ("P") found for each spectrum or containing each loss.
#'
#' @param m2l an ms2Lib Object.
#' @param to_list the type of features to search for (only P for now)
#' @param vals the type of features in which we search for the first feature (S or L for now)
#' @param reduced if TRUE, we only consider the reduced version of the patterns (with a good coverage)
#'
#' @return A list of the ids of the patterns containing the ith spectrum or loss 
#' @export
#'
#' @examples
#' data(m2l)
#' #if reduced = TRUE, calculateCoverage must have been executed first
#' list_all(m2l,"P", "S", reduced=FALSE)
#' list_all(m2l, vals = "S", reduced=FALSE)
list_all <- function(m2l, to_list = c("P"), vals=c("S", "L"), reduced = TRUE)
{
	ids <- match.arg(ids)
	if(ids == "P")
	{
		return(list_all.patterns(m2l, vals = vals, reduced = reduced))
	}
	else {
	   message("Other ids than 'P' are not implemented")
	}
}

##'@export
list_all.patterns <- function(m2l, vals = c("S", "L"), reduced=TRUE)
{
	vals <- match.arg(vals)
	if(vals == "S")
	{
		return(list_all.patterns.spectra(m2l, reduced = reduced))
	}
	else if(vals == "L")
	{
		return(list_all.patterns.losses(m2l, reduced = reduced))
	}
	else{
		message("Vals should be 'S' or 'L'")
	}
}


###select.DATAQUERIED.IDSTYPE
##' @export
list_all.patterns.spectra <- function(m2l,reduced=TRUE){
	ids <- vrange(m2l,"P",reduced=reduced)
	tp <- m2l[ids]
	sapply(patterns_from_spectra(tp,
								 length(m2l)),
		   function(x,ids){if(length(x)>0){
		   	return(ids[x])
		   	}else{
		   		return(x)
		   	}},ids=ids)
}



##' @export
list_all.patterns.losses <- function(m2l,reduced=TRUE){
	ids <- vrange(m2l,"P",reduced=reduced)
	tp <- m2l[ids]
	num_losses <- nrow(mm2EdgesLabels(m2l))
	reslist <- vector(mode = "list",length = num_losses)
	for(i in seq_along(reslist)){
		reslist[[i]] <- character(3)
	}
	idxvec <- rep(1,num_losses)

	###Replace by C++ code.
	for(i in seq_along(tp)){
		labs <- edge_attr(mm2Graph(tp[[i]]),"lab")
		for(l in labs){
			reslist[[l]] <- checkSize(reslist[[l]],idxvec[l])
			reslist[[l]][idxvec[l]] <- ids[i]
			idxvec[l] <- idxvec[l]+1
		}
	}

	###pos processing
	for(i in seq_along(reslist)){
		if(idxvec[i]>1){
			reslist[[i]] <- reslist[[i]][1:(idxvec[i]-1)]
		}else{
			reslist[[i]] <- character(0)
		}
	}

	return(reslist)
}

select.patterns.spectra <- function(m2l,id){
	temp <- select_patterns_from_spectra(mm2Patterns(m2l),as.integer(id))
	if(length(temp)>0) return(paste("P",temp,sep=""))
}


select.patterns.losses <- function(m2l,id){
	temp <- which(sapply(mm2Patterns(m2l),function(x,id){
		id %in% edge_attr(mm2Graph(x),"lab")
	},id=id))
	if(length(temp)>0) return(paste("P",temp,sep=""))
}

select.spectra.losses <- function(m2l,id){
	temp <- which(sapply(mm2Dags(m2l),function(x,id){
		id %in% edge_attr(x,"lab")
	},id=id))
	if(length(temp)>0) return(paste("S",temp,sep=""))
}


###Return f1-score,precision,recall
f1.score <- function(id,idsref,m2l,full=FALSE){

	## adds: to take the right set of spectra if the ids in the GNPS component correspond to the N indexes in the ms2Lib object 
	if("N" %in% colnames(m2l@spectraInfo))
	{
		spectres <- m2l@spectraInfo['N']
		id <- spectres[rownames(spectres) %in% id,] 
	}

	###Calculating the intersection
	inter <- intersect(id,idsref)
	if(length(inter)==0) return(rep(NA_real_,4))
	if(full&(length(inter)!=length(id))) return(rep(NA_real_,4))
	union <- union(id,idsref)

	recall <- (length(inter)/length(idsref))
	precision <- (length(inter)/length(id))
	miss_rate <- (length(idsref)-length(inter))/length(idsref)

	return(c(2*recall*precision/(recall+precision),precision,recall,miss_rate))
}




checkFTerms <- function(seq_terms){
	REF_TERMS <- c("f1","precision","recall","miss","size")
	if(all(seq_terms %in% REF_TERMS)){
		return(seq_terms)
	}else{
		stop("Unknown term(s): ",paste(seq_terms[!(seq_terms %in% REF_TERMS)],sep=", "), "authorized terms are: ",paste(REF_TERMS,sep=", "))
	}
}




#' Find the pattern with the best value of one chosen metric
#'
#' Given a list of molecules ids as characters or integers input, find the pattern maximizing one metric (F1-score, precision, recall, size) among this set of molecules.
#'
#' Should not be used directly by the user (is called by the findExplainingPatterns function)
#' 
#' @param m2l An ms2Lib object.
#' @param ids Valid IDs of objects, it mays be a character vector with the prefix "P" or an integer or numeric vector
#' @param type Which criteria is used to determine the best match, F-score, precision, recall or size of patterns (number of vertices)
#' @param returnall Shall all the patterns with a similar F1-score be returned.
#' @param full Shall only full matches be returned.
#' @param reduced Shall only the reduced set of patterns be considered.
#' @param top number of best patterns to return
#' @details The computation of the metric.
#'
#' @return A data.frame containing several slots: the ids of the elements and the values of every metric
#' @export
find.patterns.class <- function(m2l,ids,type=c("f1","precision","size"),
								returnall=FALSE,full=TRUE,reduced=FALSE, top=1){

	criterion <- checkFTerms(type)

	idsp <- vrange(m2l,"P",reduced=reduced) ## all patterns

	###The P is removed if needed.
	if(is.character(ids)&&(startsWith(ids[1],"S"))) ids <- str_sub(ids,2)

	ids <- as.numeric(ids) ## a component (clique, connected component,...) 

	###We find the best matching ones
	vf1 <- sapply(idsp,function(x,idr,m2l,fullv){
		c(f1.score(unique(m2l[x]@occurrences[,1]),idr,m2l,full=fullv),vcount(m2l[x]@graph))
	},idr=ids,m2l=m2l,fullv=full)

	rownames(vf1) <- c("f1","precision","recall","miss","size")

	###Two cases, single maximum return or the full list of ranked values returneed.
	to_return <- NULL

	###In the two cases only the non-na are considered
	pf1 <- which(!is.na(vf1["f1",]))

	if(length(pf1)==0) return(data.frame(id=character(0),f1=numeric(0),precision=numeric(0),
										 recall=numeric(0),miss=numeric(0),size=numeric(0)))
	to_return <- data.frame(id=idsp[pf1],f1=vf1["f1",pf1],
							precision=vf1["precision",pf1],recall=vf1["recall",pf1],
							miss=vf1["miss",pf1],size=vf1["size",pf1],stringsAsFactors = FALSE)
	
	to_return <- to_return[do.call(order, c(decreasing = TRUE, data.frame(to_return[,criterion]))),]
	if(!returnall){
		maxval <- to_return[1,criterion]
		posret <- 1
		while((posret < nrow(to_return)) && (!is.na(to_return[posret+1,1]))&&
			  all(to_return[posret+1,criterion] == maxval)){
			posret <- posret+1
		}
		if(posret < top)
		{
			while(posret < top && !is.na(to_return[posret+1,1]))
			{
				posret <- posret+1
			}
		}
		to_return <- to_return[1:posret,]
	}
	return(to_return)
}
