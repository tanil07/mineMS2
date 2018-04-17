


#' Default table of penalized loss
#'
#' The losses are penalized by the scoring function,
#' they are extracted form the paper by Bocker et al 2010. This includes both very
#' frequent losses, which will lead ot irllevant patterns and forbidden losses.
#'
#' @return A data frame with 2 columns the names of the losses and the
#' formula of the losses. Only the second column is used the first column is informative.
#' @export
#'
#' @examples
#' penalizedLossesDefault()
penalizedLossesDefault <- function() {
	structure(
		list(
			name = c(
				"Water",
				"Methane",
				"Ethene",
				"Ethine",
				"Butene",
				"Pentene",
				"Benzene",
				"Formaldehyde",
				"Methanol",
				"Carbon Monoxide",
				"Formic acid",
				"Carbon dioxide",
				"Acetic acid",
				"Ketene",
				"Propionic acid",
				"Malonic acid",
				"Malonic anhydride",
				"Pentose equivalent",
				"Deoxyhexose equivalent",
				"Hexose equivalent",
				"Hexuronic equivalent acid",
				"Ammonia",
				"Methylamine",
				"Methylimine",
				"Trimethylamine",
				"Cyanic Acid",
				"Urea",
				"Phosphonic acid",
				"Phosphoric acid",
				"Metaphosphoric acid",
				"Dihydrogen vinyl phosphate",
				"Hydrogen sulfide",
				"Sulfur",
				"Sulfur dioxide",
				"Sulfur trioxide",
				"Sulfuric acide",
				"Dicarbon monoxide",
				"Tetracarbon monoxide",
				"Unsaturated cyclopropane",
				"Unsaturated cyclopentane",
				"Unsaturated cycloheptane"
			),
			formula = c(
				"H2O",
				"CH4",
				"C2H4",
				"C2H2",
				"C4H8",
				"C5H8",
				"C6H6",
				"CH2O",
				"CH4O",
				"CO",
				"CH2O2",
				"CO2",
				"C2H4O2",
				"C2H2O",
				"C3H6O2",
				"C3H4O4",
				"C3H2O3",
				"C5H8O4",
				"C6H10O4",
				"C6H10O5",
				"C6H8O6",
				"NH3",
				"CH5N",
				"CH3N",
				"C3H9N",
				"CHNO",
				"CH4N2O",
				"H3PO3",
				"H3PO4",
				"HPO3",
				"C2H5O4P",
				"H2S",
				"S",
				"SO2",
				"SO3",
				"H2SO4",
				"C2O",
				"C4O",
				"C3H2",
				"C5H2",
				"C7H2"
			)
		),
		.Names = c("name", "formula"),
		class = "data.frame",
		row.names = c(NA,-36L)
	)
}




#' Discretize a set of mass differences.
#'
#' Discretize a set of mass spectra by mass losses and build a set of dags using this.
#'
#' @param m2l The ms2 lib object to be discretized.
#' @param ppm the maximum authorized deviation in ppm (parts per million).
#' @param dmz The maximum authorized deviation in Da.
#' @param freq The frequency threshold of the dataset.
#' @param limMzFormula An interval giving the range in which the formula will be calculated.
#' peaks with a masses lower than the lower term will be ignored while peak with a mass higher than the
#' higher term won't have formula check.
#' @param maxFrag The maximum number of fragment allower on the spectra.
#' @param maxOverlap The degree of overlap allowed between bins. Bins which overlap more are fused.
#' @param strictMatching Shall the formula matching be strict, or approximated.
#' @param precPpm The ppm tolerance used ot match the precursor mz to the fragments.
#' @param precDmz The minimum tolerance used to match the precursor mz to the fragments.
#' @param atoms A list of the maximum number of atoms authorized during the formula generation process.
#' If it is NULL, the default list is C:limMzFormula[2]%/%12, H:50, N:6 ,O:6.
#' @param heteroAtoms a boolean,  ignored if a custom atoms is furnished, add heteroatoms P, Cl, S as possible atoms.
#' @param mzDigits The number of decimals to write as an information in the graphs.
#' @param penalizedLosses A chracter vector, the losses which will be penalized. If NULL the table penalizedLossDefault() is used.
#'
#' @return An ms2Lib object with filled fileds and constructed dags.
#' @export
#'
#' @examples
#' print("Examples ot be put here !")
setMethod("discretizeMassLosses", "ms2Lib", function(m2l,
													 ppm = 7,
													 dmz = 0.002,
													 freq = 2/length(m2l),
													 limMzFormula = c(14.5, 200),
													 maxFrags = 15,
													 maxOverlap = 0.05,
													 strictMatching = TRUE,
													 precPpm = 20,
													 precDmz = 0.02,
													 atoms = NULL,
													 heteroAtoms = TRUE,
													 mzDigits = 4,
													 penalizedLosses = NULL,
													 majoredLosses = character(0),
													 ...) {

	message("Discretization of the mass losses...")

	###Paameters checking
	if(maxFrag>20){
		stop("mineMS2 is not made to process more than 20 fragments by spectra, please lower
			 the maxFrag parameters")
	}
	if(maxOverlap>0.15){
		warning("maxOverlap is superior to 0.15, this can lead to incoherent labels.")
	}

	if((!is.numeric(limMzFormula))|(length(limMzFormula)!=2)){
		stop("limMzFormula should be a range with 2 floats")
	}

	if((limMzFormula[1]<=10|limMzFormula[2]>250)){
		limMzFormula <- c(max(limMzFormula[1],10),min(max(limMzFormula[2],250)))
		warning("limMzFormula is to wide, to avoid computation issue it have been set to ",
				 limMzFormula[1],"-",limMzFormula[2])
	}


	if(is.null(atoms)){
		if(heteroAtoms){
			atoms <- list("C"= max(limMzFormula)%/%12,"H"=50,"N"=6,"O"=6)
		}else{
			atoms <- list("C"=max(limMzFormula)%/%12,"H"=50,"N"=6,"O"=6,"S"=2,"Cl"=1,"P"=2)
		}
	}


	if(is.null(penalizedLosses)){
		penalizedLosses <- penalizedLossesDefault()$formula
		message("penalizedLosses not furnished, default table is used.")
	}else{
		penalizedLosses <- sapply(penalizedLosses,stringToFormula,vnames=names(l_atoms),simplify=FALSE)
		penalizedLosses <- sapply(penalizedLosses,formulaToString)
	}

	res_list <- discretizeMassesDifferences(mm2Spectra(m2l),
								ppm = ppm, dmz = dmz,
								freq = freq, mzdigits = mzDigits,
								limFormula = limMzFormula, maxFrag = maxFrag,
								max_overlap = maxOverlap, strictMatching = strictMatching,
								prec.ppm = precPpm,prec.dmz = precDmz, atoms = atoms,
								floss = penalizedLosses, nfloss = majoredLosses, ...)

	###Add the fusing part of the edge labels.
	message("Fusing the edges labels.")

	###Simplifying the DAG if necessary.
	change <- TRUE
	niter <- 1
	while(change){
		res_list <- fuseElem(elems =res_list$elems ,dags = res_list$dags,atoms = names(atoms))
		change <- res_list$change
		niter <- niter+1
	}


	mm2EdgesLabels(m2l) <- res_list$elem
	mm2Dags(m2l) <- res_list$dags
	message("Discretization finished.")
	m2l
})


###Return an extended spectra and possibly add the precurosr if needed.
###If maxnum is not NULL return the top "k" or top "k-1" depending of the inclusion
###of the precursor.
#' Return the spectrum considered while building the graph.
#'
#' @param spec A "Spectrum2" object as defined in msnBase
#' @param ppm A ppm parameter to match the precursor
#' @param dmz A minimum deviation in dalton to match the precursor.
#' @param removeSup A boolean indicating if peaks with a mass higher than the precursor should be remode
#' @param addPrecursor A boolen indicating if a preucrsor should be added if needed
#' @param maxNum The maximum number of peak to keep. The value correspond ot the top-k intensities.
#' @param relInt Should the relative intensity be added to the dataset.
#'
#' @return A data.frame containg the mz, intensity and precursor columns.
#' @export
#'
#' @examples
#' print('examples to be put there')
getExtendedSpec <- function(spec,ppm=5,dmz=0.02,removeSup=TRUE,addPrecursor=TRUE, maxNum = NULL, relInt = TRUE){



	###Reconstruct a spectra form an msnbase object.
	df_spec <- data.frame(mz=mz(spec),intensity=intensity(spec))

	if(nrow(df_spec)==0){
		return(data.frame(mz=numeric(0),intensity=numeric(0),precursor=logical(0)))
	}

	mzprec <- precursorMz(spec)
	tol <- max(mzprec*ppm*(10e-6),dmz)

	###We remove the values with a masses higher than the precursor.
	if(removeSup){
		psup <- which(df_spec$mz<(mzprec+tol))
		if(length(psup)>0) df_spec <- df_spec[psup,,drop=FALSE]
	}

	###intensity is the minimal intensity of the spectra if the peak does not exists.
	intprec <- min(df_spec$int)

	###Checking if there is a potential precursor
	vec_prec <- NULL
	pprec <- which(abs(df_spec$mz-mzprec)<tol)
	if(length(pprec)==0){
		vec_prec <-rep(FALSE,nrow(df_spec))
		if(addPrecursor){
			df_spec[nrow(df_spec)+1,] <- c(mzprec,intprec)
			vec_prec[nrow(df_spec)] <- TRUE
		}
	}else{
		vec_prec <- rep(FALSE,nrow(df_spec))
		vec_prec[pprec[which.min(df_spec$mz[pprec]-mzprec)]]<-TRUE
	}
	df_spec$precursor <- vec_prec
	pos_prec <- which(vec_prec)

	###Precursor is included in every case.
	if(!is.null(maxNum)){
		if(maxNum < nrow(df_spec)){
			pm <- order(df_spec$int,decreasing=TRUE)
			tpos <- pm[1:maxNum]
			if(addPrecursor & !(pos_prec %in% tpos)){
				tpos <- c(tpos[1:(maxNum-1)],pos_prec)
			}
			df_spec <- df_spec[tpos,]
		}
	}

	if(relInt){
		onames <- colnames(df_spec)
		onames <- c(onames,"rel_int")
		rel_int <- (as.numeric(df_spec$int)/sum(as.numeric(df_spec$int)))*100
		if(any(is.na(rel_int))) browser()
		df_spec$rel_int <- rel_int
	}
	df_spec <- df_spec[order(df_spec$mz,decreasing=TRUE),]
	df_spec
}

###Discretize set of values furnished in the vmz vector
###with a bandwidth given as the minimam between ppm and dmz.
###frac give the minimum number of samples.
discretizeSequenceByClasses <- function(ppm,
			 vmz, vsample, ndiff, frac = 0,
			 inter = 0.2, dmz = 0.01,
			 nPoints = 512, extended=FALSE) {
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
			plim = findLimDensity(den$y, plim[2] + 1, plim[3])
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

##Generate a matrix of differences given a sequence.
generateMatDiff <- function(seqmz){
	rmat <- outer(seqmz, seqmz, FUN = "-")
	rmat[lower.tri(rmat, diag = TRUE)] <- NA
	return(list(
		mz = rmat[upper.tri(rmat)],
		idx = which(upper.tri(rmat), arr.ind = TRUE)
	))
}


convertMatrixIgraph <- function(x){
	x[which(!is.na(x),arr.ind = TRUE)] <- 1
	x[which(is.na(x),arr.ind = TRUE)] <- 0
	x
}


checkFracParam <- function (num, maxVal, OneEqualAll = TRUE)
{
	if (num != round(num)) {
		if (num > 1) {
			stop("Param too high.")
		}
		else {
			return(floor(maxVal * num))
		}
	}
	else {
		if (num > 1 & num <= maxVal)
			return(num)
		if (num > maxVal)
			stop("Param too high.")
		if (num == 1 & OneEqualAll)
			return(maxVal)
		if (num == 1)
			return(1)
		if (num == 0)
			return(0)
	}
}


sumFormula <- function(f1,f2,all_lab=NULL){
	if(is.na(f1)|is.na(f2)) return(NA)
	if(is.null(all_lab)){
		all_lab <- MS2process:::tabAtoms()$name
	}
	vf <- rep(0,length(all_lab))
	names(vf) <- all_lab
	vf[all_lab] <- vf[all_lab]+f1
	vf[all_lab] <- vf[all_lab]+f2
	###We remover the 0
	vf <- vf[vf!=0]
	return(vf)
}



discretizeMassesDifferences <- function(list_spec,
										ppm = 7, dmz = 0.002,
										freq = 0.05, limFormula=c(14.5,200),
										mzdigits=3, maxFrag = 25,
										max_overlap=0.05,strictMatching=TRUE,
										prec.ppm=20,prec.dmz=0.02,
										atoms = list("C"= max(limMzFormula)%/%12,"H"=50,"N"=6,"O"=6) ,
										floss=NULL,nfloss = NULL,...) {
	list_spec <- lapply(list_spec,getExtendedSpec,maxNum=maxFrag,ppm=prec.ppm,
						dmz=prec.dmz,relInt = TRUE)
	list_masses <- lapply(list_spec,'[[',i="mz")
	list_int <- lapply(list_spec,'[[',i="rel_int")

	mean_mzf <- sum(limFormula)/2
	tol <- diff(limFormula)/2
	l_atoms <- atoms

	if(is.null(floss)){
		data("freq_loss")
		floss <- sapply(freq_loss[,2],stringToFormula,vnames=names(l_atoms),simplify=FALSE)
		floss <- sapply(floss,formulaToString)
	}
	if(is.null(nfloss)){
		data("non_freq_loss")
		nfloss <- sapply(non_freq_loss[,2],stringToFormula,vnames=names(l_atoms),simplify=FALSE)
		nfloss <- sapply(nfloss,formulaToString)
	}

	ndiff <- length(list_masses)
	thresh <- checkFracParam(freq, ndiff)
	list_matrix <- lapply(list_masses,generateMatDiff)
	pcord <- lapply(list_matrix, '[[', i = "idx")
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

	###Now we create all the standard deviation of the masses.
	masses <- resdisc$discrete_elem[,1]
	diffv <- resdisc$discrete_elem[,3]-resdisc$discrete_elem[,2]

	###We define the max larg
	qv <- resdisc$discrete_elem[,4]

	###The standard deviation is created
	sds <- qv

	fac_sig <- 2
	###Merging overlapping values.
	message("Masses merging\n")

	merged_masses <- gaussianMerging(masses,sds^2,alpha=max_overlap,fac_sig = fac_sig)
	high_mz_idx <- which(merged_masses$mu>limFormula[2])
	###Now we check the values which may originates from a formula.


	###Generation of all the neutral formula
	message("Formula generation\n")

	###How ot handle heteroatoms, maximum nuimber of different heteroatoms.
	allFormula <- lossesFormulaGeneration(limFormula,atoms=l_atoms,...)

	##Associating formula with losses.
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
			dag <- add_edges(dag,as.numeric(t(pf)),lab=toReturn[[i]][pf])
			###We increment the counter value.
			for(inc in toReturn[[i]][pf])	recap_tab$count[inc] <- recap_tab$count[inc] +1
		}
		toReturn[[i]] <- dag
	}
	return(list(spectra=list_spec,dags = toReturn, elems = recap_tab))
}



###Function used for labels fusion if there is an incoherence.
fuseElem <- function(elems,dags,atoms=NULL){
	if(is.null(atoms)) atoms <- tabAtoms()$name
	anyChange <- FALSE

	###We first store a list of all the formula when they exist
	allformula <- str_split(elems$formula,fixed("|"))
	allformula <- sapply(allformula,function(x,vi){
		if(is.na(x[1])|x[1]=="NA") return(NULL)
		sapply(x,function(y,vn){
			stringToFormula(y,vnames = vn)
		},vn=vi,simplify=FALSE)
	},vi=atoms,simplify=FALSE)


	###Data structure for triple value, a sparse matrix.
	stormat <- Matrix(0,nrow=nrow(elems),ncol=nrow(elems),sparse = TRUE)
	storval <- vector(mode="list",length=1000)
	nval <- 1

	####We get a list of all the possible error in the labels.
	for(igl in 1:length(dags)){
		g <- dags[[igl]]
		if(is.null(g)||is.na(g)){
			next
		}
		adj_list <- adjacent_vertices(g,V(g),"out")
		lab_sec <- E(g)$lab
		for(n in V(g)){
			adjv <- adj_list[[n]]
			if(length(adjv)>0){
				for(na in adjv){
					succ <- intersect(adj_list[[n]],adj_list[[na]])
					if(length(succ)!=0){
						for(s in succ){
							all_ids <- get.edge.ids(g,c(n,na,na,s,n,s))
							la <- lab_sec[all_ids[1]]
							lb <- lab_sec[all_ids[2]]
							lc <- lab_sec[all_ids[3]]

							ma <- min(la,lb)
							mb <- max(la,lb)

							###We check if the value may be added ot the dataset or not.
							is_in <- stormat[ma,mb]
							if(is_in==0){
								###If necessary the storage vector is pushed up
								if(nval>length(storval)){
									storval <- c(storval,vector(mode="list",length=length(storval)))
								}
								###The position is added ot the dataset.
								stormat[ma,mb] <- nval
								storval[[nval]] <- c(la,lb,lc)
								nval <- nval+1
							}else{
								if(!lc %in% storval[[is_in]] ){
									storval[[is_in]] <- c(storval[[is_in]],lc)
								}
							}
						}
					}
				}
			}
		}
	}
	storval <- storval[1:(nval-1)]

	####Now we have a data structure with all the value, we check which set needs to be fused.
	vmul <- sapply(storval,function(x){length(unique(x))>3})
	pmistake <- which(vmul)

	nreduc <- 0
	if(length(pmistake)>0){
		anyChange <- TRUE
		nreduc <- sapply(storval[pmistake],length)-3
	}

	warning("There is ",length(vmul),"incompatibilities, this may be caused by a too high accuraccy in mass.",sum(nreduc),"labels will be removed.")

	###We we fuse the dataset if necessary
	to_rm <- numeric(0)
	nlab <- elems$lab
	merged_val <- numeric(0)

	for(i in pmistake){
		sv <- storval[[i]]

		a <- sv[1]
		b <- sv[2]
		cs <- sv[3:length(storval[[i]])]
		telems <- elems[cs,,drop=FALSE]

		##The dataset are fused.
		mzmin <- min(telems$mzmin)
		mzmax <- max(telems$mzmax)
		mz <- mean(telems$mz)
		fused <- TRUE
		count <- sum(telems$count)

		###We determine correct formula.
		tform <- allformula[cs]

		formula <- NULL

		non_null <- which(!sapply(tform,is.null))
		if(length(non_null)==0){
			formula <- "NA"
		}else{
			###We calculate the new formula.
			aformula <- ifelse(is.null(allformula[[a]]),NA,allformula[[a]])
			if(typeof(aformula)=="list"){
				aformula <- aformula[[1]]
			}
			bformula <- ifelse(is.null(allformula[[b]]),NA,allformula[[b]])
			if(typeof(bformula)=="list"){
				bformula <- bformula[[1]]
			}
			cformula <- sumFormula(aformula,bformula,all_lab = atoms)

			if(!is.na(cformula)){
				resf <- formulaToString(cformula)

				cform <- unlist(tform,recursive = FALSE,use.names = TRUE)

				ncform <- unique(names(cform))
				pfor <- match(resf,ncform)

				if(is.na(pfor)){
					formula <- c(resf,ncform)
				}else{
					formula <- c(resf,ncform[-pfor])
				}
			}else{
				formula <- "NA"
			}
		}
		formula <- paste(formula,sep="|")
		to_rm <- c(to_rm,cs[2:length(cs)])
		merged_val <- c(merged_val,cs[1])
		elems$mzmin[cs[1]] <- mzmin
		elems$mzmax[cs[1]] <- mzmax
		elems$mz[cs[1]] <- mz
		elems$fused[cs[1]] <- TRUE
		elems$count[cs[1]] <- count
		elems$formula[cs[1]] <- formula
		nlab[cs] <- cs[1]
	}

	###Now all the labels need to be redone
	ulabs <- unique(nlab)
	new_labs <- match(nlab,ulabs)

	###Relabeling of the graph
	for(igl in 1:length(dags)){
		g <- dags[[igl]]
		if(is.null(gl)||is.na(gl)){
			next
		}
		eg_lab <- edge_attr(g,"lab")
		edge_attr(g,"lab") <- new_labs[eg_lab]
		dags[[igl]] <- g
	}

	###Updating the lement table.
	to_rm <- to_rm[!(to_rm %in% merged_val)]


	if(length(to_rm)>=1){
		elems <- elems[-to_rm,]
	}
	elems$lab <- 1:nrow(elems)
	return(list(elem=elems,dag=dags,change=anyChange))
}
