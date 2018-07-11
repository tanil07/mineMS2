####Implementation of the selecte functions.


checkSize <- function(vec,pos){
	if(length(vec)<pos){
		vec <- c(vec,do.call(typeof(vec),args=list(length(vec))))
	}
	return(vec)
}

#' Selection function
#'
#' @param m2l An ms2Lib object.
#' @param ids Valid IDs of objects.
#'
#' @return The objects of the correspodning type.
#' @export
#'
#' @examples
#' print("Examples to be put here")
select <- function(m2l,ids,vals=c("P","S","L","F"),...){
	if(class(m2l)!="ms2Lib") stop("m2l should be an ms2Lib object.")

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

#' find patterns or spectra form a set of pattern or spectra ids
#'
#' find data of type "vals" containing all the ids listed in "ids"
#'
#' @param m2l an ms2Lib Object.
#' @param ids A list of ids to be searched for.
#' @param vals The type of the returned values.
#' @param ... supp arguments.
#'
#' @return The ids of the patterns or spectra cotaining all the ids.
#' @export
#'
#' @examples
#' print("examples to be put here")
find <- function(m2l,ids,vals=c("P","S"),...){
	if(class(m2l)!="ms2Lib") stop("m2l should be an ms2Lib object.")

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


find.patterns.spectra<- function(m2l,ids,reduced=FALSE,biggest=TRUE,partial = TRUE){
	idsp <- vrange(m2l,"P",reduced)

	resids <- character(10)
	resids_partial <- character(10)
	seqsize <- integer(10)
	seqsize_partial <- integer(10)
	nsol <- 0
	nsol_partial <- 0
	for(ip in seq_along(idsp)){
		temp_pat <- m2l[idsp[ip]]
		occs <- mm2Occurences(temp_pat)

		if(all(ids %in% occs[,"gid"])){
			nsol <- nsol+1
			resids <- checkSize(resids,nsol)
			resids[nsol] <- idsp[ip]
			seqsize <- checkSize(seqsize,nsol)
			seqsize[nsol] <- vcount(temp_pat@graph)-1

		}

		if(any(ids %in% occs[,"gid"])){
			nsol_partial <- nsol_partial+1
			resids_partial <- checkSize(resids_partial,nsol_partial)
			resids_partial[nsol_partial] <- idsp[ip]
			seqsize_partial <- checkSize(seqsize_partial,nsol_partial)
			seqsize_partial[nsol_partial] <- sum(ids %in% occs[,"gid"])+vcount(temp_pat@graph)-1
		}
	}

	if(nsol!=0){
		resids <- resids[1:nsol]
		if(biggest){
			resids <- resids[which.max(seqsize)]
		}
		return(resids)
	}else{
		if(partial){
			if(length(resids_partial)>0){
				resids_partial <- resids_partial[which.max(seqsize_partial)]
			}
		}
		return(character(0))
	}
}



###select.DATAQUERIED.IDSTYPE

all.patterns.spectra <- function(m2l,reduced=TRUE){
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




all.patterns.losses <- function(m2l,reduced=TRUE){
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
	temp <- select_patterns_from_spectra(mm2Patterns(m2l),as.integer(substring(id,2)))
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

