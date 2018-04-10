
#' Discretization of mass differences.
#'
#' Creates a discretized set of discrtized edges labels.
#'
#' @export
#' @param x poney
#' @param supp_infos Supplementary information to be associated to the spectra.
#' @return no valus is returned. To access the set of labels use...
#' @examples
#' print("examples to be put here")

setMethod("discretizeEdges","ms2Lib", function(m2l){

})



###Discrtize a set of mass differneces using a gaussian density.
discretizeSequenceByClasses <-
	function(ppm,
			 vmz,
			 vsample,
			 ndiff,
			 frac = 0,
			 inter = 0.2,
			 dmz = 0.01,
			 nPoints = 512,
			 rlgroup = TRUE,extended=FALSE) {
		porder <- order(vmz)
		vall <- vmz[porder]
		mzrange <- range(vall)
		mzrange <-
			c(max(floor(mzrange[1] / 10) * 10, 1), ceiling(mzrange[2] / 10) * 10)
		massInter <- seq(mzrange[1], mzrange[2], inter / 2)
		masspos <- FindEqualGreaterM(vall, massInter)
		num_group <- 0
		idxgroup <- vector(mode = "list", 512)
		listgroup <- matrix(0, nrow = 512, ncol = 4)
		pos = 0
		previousMz = NULL
		num_group = 0
		for (i in seq(length = length(massInter) - 2)) {
			mem = i
			start <- masspos[i]
			end <- masspos[i + 2] - 1
			if (end - start < 0)
				next
			subpeakl = vall
			pos = pos + 1
			bw = max((ppm * massInter[i] * 1e-6), dmz)/3
			den <-
				density(
					subpeakl,
					bw,
					from = massInter[i] - 3 * bw,
					to = massInter[i + 2] +
						3 * bw,
					n = nPoints
				)

			###Putting the close to 0 to 0
			maxy = max(den$y)
			den$y[which(den$y <= bw * maxy)] = 0
			plim = c(-1, 0, 0)
			if (den$y[2] > den$y[1])
				plim[3] <- 1
			oldMz = previousMz
			repeat {
				plim = proFIA:::findLimDensity(den$y, plim[2] + 1, plim[3])
				if (plim[2] <= plim[1])
					break
				selectedPGroup = which(subpeakl >= den$x[plim[1]] &
									   	subpeakl <= den$x[plim[2]])
				if (length(selectedPGroup) == 0)
					next
				num_group = num_group + 1
				if (num_group > nrow(listgroup)) {
					listgroup <- rbind(listgroup,
									   matrix(
									   	0,
									   	nrow = nrow(listgroup),
									   	ncol = 4
									   ))
				}
				if (mean(vmz[porder[selectedPGroup]]) %in% oldMz ||
					(abs(median(subpeakl[selectedPGroup]) - massInter[i]) < bw *
					 2) ||
					(abs(median(subpeakl[selectedPGroup]) - massInter[i + 2]) <
					 bw * 2)) {
					num_group = num_group - 1
					next
				}
				if (length(unique(vsample[porder[selectedPGroup]])) / ndiff < frac) {
					num_group = num_group - 1
					next
				}
				if(extended){
					listgroup[num_group, ] <-
						c(mean(vmz[porder[selectedPGroup]]),
						  den$x[plim[1]],
						  den$x[plim[2]],bw)
				}else{
					listgroup[num_group, ] <-
						c(mean(vmz[porder[selectedPGroup]]),
						  range(vmz[porder[selectedPGroup]]),bw)
				}
				idxgroup[[num_group]] <- porder[selectedPGroup]
				previousMz <- c(previousMz, mean(vmz[porder[selectedPGroup]]))
			}
		}

		if (num_group < 1) {
			return(list(dicrete_elem = numeric(0), idx = list()))
		}
		listgroup <- listgroup[1:num_group,]
		idxgroup <- idxgroup[1:num_group]
		return(list(discrete_elem = listgroup, idx = idxgroup))
	}




convertMatrixIgraph <- function(x){
	x[which(!is.na(x),arr.ind = TRUE)] <- 1
	x[which(is.na(x),arr.ind = TRUE)] <- 0
	x
}


formulaToString <- function(vformula,vnames = NULL){
	if(is.null(vnames))	vnames <- names(vformula)
	if(is.factor(vnames)) vnames <- as.character(vnames)
	vnums <- ifelse(vformula>1,vformula,"")
	vnames <- ifelse(as.character(vformula)=="0","",vnames)
	paste(paste(vnames,vnums,sep=""),collapse="")
}


###
#1) discretize the mzs
#2) check the signals which are too close and put everyone
#3) get all the

###Test of the number of generated formula.
library(MS2process)
discretizeMassesDifferences <- function(list_spec,
										ppm = 7,
										dmz = 0.002,
										freq = 0.05,
										mzdigits=3,intdigits=2,
										limFormula=c(14.5,200), maxFrag = 25,
										max_overlap=0.02,strictMatching=TRUE,
										tolPrec=0.008, heteroAtoms=TRUE,floss=NULL,
										nfloss = NULL,...) {
	list_spec <- lapply(list_spec,getExtendedSpec,maxNum=maxFrag,tol=tolPrec)
	list_masses <- lapply(list_spec,'[[',i="mz")
	list_int <- lapply(list_spec,'[[',i="rel_int")

	mean_mzf <- sum(limFormula)/2
	tol <- diff(limFormula)/2
	l_atoms <- NULL
	if(heteroAtoms){
		l_atoms <- list(C=floor((mean_mzf+tol)/12),H=50,N=6,O=6,S=2,Cl=1,P=2)
	}else{
		l_atoms <- list(C=floor((mean_mzf+tol)/12),H=50,N=6,O=6)
	}

	if(is.null(floss)){
		data("freq_loss")
		floss <- sapply(freq_loss[,2],MS2process:::stringToFormula,vnames=names(l_atoms),simplify=FALSE)
		floss <- sapply(floss,MS2process:::formulaToString)
	}
	if(is.null(nfloss)){
		data("non_freq_loss")
		nfloss <- sapply(non_freq_loss[,2],MS2process:::stringToFormula,vnames=names(l_atoms),simplify=FALSE)
		nfloss <- sapply(nfloss,MS2process:::formulaToString)
	}

	ndiff <- length(list_masses)
	thresh <- MS2process:::checkFracParam(freq, ndiff)
	list_matrix <- lapply(list_masses,generateMatDiff)
	pcord <- lapply(list_matrix, '[[', i = "idx")
	#browser()
	vecsample <- rep(seq_along(pcord), times = sapply(pcord, nrow))
	vecsize <- sapply(pcord, nrow)
	pcord <- do.call('rbind', pcord)
	vmz <- unlist(sapply(list_matrix, '[[', i = "mz"))
	allEdges<- NULL
	valinter <- c(0, cumsum(vecsize))
	treatInter <- function(x) {
		pvec <- .bincode(x, valinter, right = TRUE)
		if (length(pvec) == 1) {
			tr <- matrix(c(pvec, pcord[x, ]), nrow = 1)
			colnames(tr) <- c("pvec", "row", "col")
			return(tr)
		}
		cbind(pvec, pcord[x, ])
	}

	cat("First discretization\n")
	###We start by discretizing all the values
	resdisc <-
		discretizeSequenceByClasses(
			ppm = ppm,
			vmz = vmz,
			vsample = vecsample,
			dmz = dmz,
			ndiff = length(list_masses),
			frac = freq
		)

	# df1 <- cbind(resdisc$discrete_elem)
	###Now we create all the standard deviation of the masses.
	masses <- resdisc$discrete_elem[,1]
	diffv <- resdisc$discrete_elem[,3]-resdisc$discrete_elem[,2]
	###We define the max larg
	qv <- resdisc$discrete_elem[,4]
	###The standard deviation is created
	sds <- qv
	# sds <- ifelse(tval>dmz,tval,dmz)/2

	fac_sig <- 2
	###Merging overlapping values.
	cat("Masses merging\n")
	merged_masses <- gaussianMerging(masses,sds^2,alpha=max_overlap,fac_sig = fac_sig)
	# df2 <- cbind(merged_masses)
	high_mz_idx <- which(merged_masses$mu>limFormula[2])
	###Now we check the values which may originates from a formula.


	###Generation of all the neutral formula.



	cat("Formula generation\n")
	# allNeutral <- findFormulaCpp(mean_mzf,tol,rules=c(FALSE,FALSE,rep(FALSE,5)),
	# 							 atoms = l_atoms,mgraph=FALSE)

	###How ot handle heteroatoms, maximum nuimber of different heteroatoms.
	allFormula <- RformulaExtension(limFormula,atoms=l_atoms,...)




	###We calculate all the masses.
	# df_at <- MS2process:::makeDfValidAtoms(l_atoms,limFormula[2])
	# f_masses <- sapply(allNeutral,function(x,tatoms){
	# 	sum(tatoms$mass*x)
	# },tatoms=df_at)

	##New version
	f_masses <- allFormula$masses

	###Now we generate the bmin to fac_sig
	tval <- f_masses*ppm/10^6
	sd_f <- ifelse(tval>dmz,tval,dmz)

	maxm <- f_masses+fac_sig*sd_f
	minm <- f_masses-fac_sig*sd_f

	merginf <- merged_masses$mu-fac_sig*sqrt(merged_masses$sig)
	mergmax <- merged_masses$mu+fac_sig*sqrt(merged_masses$sig)

	###Now we need to check the overlap between the interval.
	##Function to parse the Cpp function.
	cat("Interval checking\n")
	# resinter <- Rcpp::sourceCpp('~/ms-ms/Cpp/IntervalMatching.cpp')
	matched_inter <- checkInter(minm,maxm,merginf,mergmax)
	vec_inter <- do.call("cbind",matched_inter)

	atm_names <- colnames(allFormula$formula)
	###Now we create a table centralizing all the values.
	str_formula <- NULL
	freq_f <- NULL
	non_freq_f <- NULL

	xp1 <- 1
	xp2 <- 2
	if(strictMatching){
		xp1 <- 3
		xp2 <- 4
	}

	str_formula <- apply(vec_inter,1,function(x){
		if(is.na(x[xp1])) return(NA_character_)
		apply(allFormula$formula[x[xp1]:x[xp2],,drop=FALSE],1,formulaToString,vnames=atm_names)
	})
	###Finding penalized and advantaged formula.
	freq_f <- sapply(str_formula,function(x){
		if(is.na(x)) return(FALSE)
		any(x %in% floss)
	})

	non_freq_f <- sapply(str_formula,function(x){
		if(is.na(x)) return(FALSE)
		any(x %in% nfloss)
	})


	###Final strings
	str_formula <- sapply(str_formula,function(x){
		if(is.na(x)) return(NA_character_)
		paste(x,collapse = "|")
	})

	####We check if the losses are nitrogen or carbon only.
	carb_only <- str_detect(str_formula,"^C[0-9]{0,10}$")
	nitrogen_only <- str_detect(str_formula,"^N[0-9]{0,10}$")
	CH_only <- str_detect(str_formula,"^C[0-9]{0,10}H[0-9]{0,10}$")

	###we filter out intervals without formula.
	to_keep <- which(!is.na(str_formula))

	###Now we build a recapitulative table, keeping only the mean mz,min mz, formula and label.
	recap_tab <- data.frame(lab=1:(length(to_keep)+length(high_mz_idx)),mz=merged_masses$mu[c(to_keep,high_mz_idx)],
							mzmin=merginf[c(to_keep,high_mz_idx)],
							mzmax=mergmax[c(to_keep,high_mz_idx)],
							sig=merged_masses$sig[c(to_keep,high_mz_idx)],
							fused = merged_masses$fused[c(to_keep,high_mz_idx)],
							count=numeric((length(to_keep)+length(high_mz_idx))),
							formula=c(str_formula[to_keep],rep(NA_character_,length(high_mz_idx))),
							adv_loss = c(freq_f[to_keep],rep(FALSE,length(high_mz_idx))),
							pen_loss = c(non_freq_f[to_keep],rep(FALSE,length(high_mz_idx))),
							carb_only = c(carb_only[to_keep],rep(FALSE,length(high_mz_idx))),
							nitrogen_only = c(nitrogen_only[to_keep],rep(FALSE,length(high_mz_idx))),
							ch_only = c(CH_only[to_keep],rep(FALSE,length(high_mz_idx))),
							stringsAsFactors = FALSE)

	# browser()

	###Graph creations
	toReturn <- vector(mode = "list", length = length(list_masses))




	for (i in seq_along(list_masses)){
		toReturn[[i]] <- matrix(
			NA_real_,
			nrow = length(list_masses[[i]]),
			ncol = length(list_masses[[i]])
		)
		colnames(toReturn[[i]]) <- sprintf(paste('%0.',mzdigits,'f',sep=""),list_masses[[i]])

		###We first need to reorder the mass diff
		om_diff <- order(list_matrix[[i]]$mz,decreasing = FALSE)
		list_matrix[[i]]$mz <- list_matrix[[i]]$mz[om_diff]
		list_matrix[[i]]$idx <- list_matrix[[i]]$idx[om_diff,,drop=FALSE]

		matched <- disjointBins(list_matrix[[i]]$mz, recap_tab$mzmin, recap_tab$mzmax, recap_tab$mz)
		###We keep only the matched edges.
		pok <- which(!is.na(matched))

		if(length(pok)==0){
			toReturn[[i]]<-NA
			next
		}
		###We keep only the labeled edges.
		list_matrix[[i]]$mz <- list_matrix[[i]]$mz[pok]
		list_matrix[[i]]$idx <- list_matrix[[i]]$idx[pok,,drop=FALSE]
		toReturn[[i]][list_matrix[[i]]$idx]<-matched[pok]
		###We create the graph.
		idxnodes <- match(list_masses[[i]],sort(list_masses[[i]],decreasing = FALSE))-1
		dag <- graph.empty(n = nrow(toReturn[[i]]),directed=TRUE)
		dag <- set_vertex_attr(dag,name="mz",value = list_masses[[i]])
		dag <- set_vertex_attr(dag,name="rel_int",value = list_int[[i]])
		dag <- set_vertex_attr(dag,name="label",value = as.integer(idxnodes))
		dag <- set_vertex_attr(dag,name="prec",value=as.character(list_spec[[i]]$prec))
		pf <- which(!is.na(toReturn[[i]]),arr.ind = TRUE)
		#pf <- pf[,c(2,1)]
		to_remove <- numeric(100)
		nrm <- 1
		if(nrow(pf)!=0){
			##We remove the parent with the same descendant.
			lab_entry <- paste(pf[,2],toReturn[[i]][pf],sep="_")

			pdup <- table(lab_entry)

			to_check <- which(as.numeric(pdup)>=2)
			if(length(to_check)!=0){
				for(j in to_check){
					if(nrm>length(to_remove)) ro_remove <- c(to_remove,numeric(length(to_remove)))
					###We get the concerned edges
					posf <- which(lab_entry==names(pdup)[j])

					prob_edges <- pf[posf,]
					###Now we calculate the masses.
					dmz <- list_masses[[i]][prob_edges[,1]]-list_masses[[i]][prob_edges[,2]]
					mass_to_comp <- recap_tab[toReturn[[i]][prob_edges][1],"mz"]
					if(all(is.na(mass_to_comp))) next

					pmin <- which.min(abs(dmz-mass_to_comp))
					if(length(pmin)!=0){
						to_remove[nrm:(nrm+length(posf)-2)] <- posf[-pmin]
					}
					nrm <- nrm + length(posf)-1
				}
			}

			to_remove <- to_remove[1:(nrm-1)]
			if(nrm>1){
				pf <- pf[-to_remove,]
			}
			# print(pf)
			dag <- add_edges(dag,as.numeric(t(pf)),lab=toReturn[[i]][pf])
			###We increment the counter value.
			for(inc in toReturn[[i]][pf])	recap_tab$count[inc] <- recap_tab$count[inc] +1
		}
		toReturn[[i]] <- msDag(dag,list_spec[[i]])
	}
	return(list(dag=toReturn,elem = recap_tab))
}
