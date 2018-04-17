

#' findFormula function
#'
#'Fast formula generator using the algoirthm descirbed by Bocker
#'2010. 100 times faster than the mgf used.
#' @param mz The mass to b decomposed.
#' @param tol The tolerance in ppm.
#' @param atoms A of constraint on the number of atoms.
#'
#' @return A list of the valid formula.
#' @export
#'
#' @examples
#' cat("examples ot be put there")
findFormula <-
	function(mz,
			 tol,
			 atoms = list(H = 50,
			 			 C = 25,
			 			 N = 10,
			 			 O = 20),
			 rules = rep(TRUE, 7),
			 mgraph = TRUE,
			 scomp) {
		df <- makeDfValidAtoms(atoms, mz)
		b <- 200
		allSol <- decomposeMass(mz,tol,as.integer(b),df$mass)
		to_keep <- numeric(0)
		if (length(allSol) == 0) {
			return(list())
		}
		##We filter out any solution which does not satisfy the criterion
		to_rm <- rep(FALSE, length(allSol))
		for (i in seq_along(allSol)) {
			for (j in 1:length(atoms)) {
				if (allSol[[i]][j] > atoms[[j]]) {
					to_rm[i] <- TRUE
					break
				}
			}
		}
		to_keep <- which(!to_rm)
		if (length(to_keep) == 0)
			return(list())
		allSol <- allSol[to_keep]
		to_keep <-
			sapply(
				allSol,
				checkSevenGoldenRules,
				tAtoms = df,
				tol = 'high',
				rules = rules
			)
		if (length(to_keep) > 0) {
			allSol <- allSol[to_keep]
		}
		if (length(allSol) == 0) {
			return(list())
		}
		resList <- vector(mode = "list", length = length(allSol))
		nsol <- 1
		for (i in 1:length(allSol)) {
			###We check that the formula should be kept.
			rb <- checkGraph(allSol[[i]],
							 t_atoms = df,
							 singleComponent = scomp)
			if (rb | !(mgraph)) {
				resList[[nsol]] <- allSol[[i]]
				names(resList[[nsol]]) <- df[, "name"]
				nsol <- nsol + 1
			}
		}
		return(resList[1:nsol])
	}



###Reference for atoms mass etc...
tabAtoms <- function(){
	n_atoms <- c("H", "C", "O", "N", "F", "Cl", "Br", "I", "Si",
				 "S", "P", "Se")
	m_atoms <- c(1.00783, 12, 15.994915, 14.003074, 18.998403,
				 34.0968853, 78.918336, 126.904477, 27.976928, 31.972071,
				 30.973763, 73.922477)
	val_atoms <- c(1, 4, 2, 3, 1, 1, 1, 1, 4, 2, 3, 6)
	return(data.frame(name = n_atoms, mass = m_atoms, valence = val_atoms,
					  stringsAsFactors = FALSE))
}

###Construct a reference table form a lsit of atoms.
makeDfValidAtoms <- function (clist, mz = NULL)
{
	tAtoms <- tabAtoms()
	ref <- names(clist)
	pok <- ref %in% tAtoms$name
	if (any(!pok)) {
		stop("invalid atoms furnished :", paste(ref[!pok], collapse = ",",
												sep = " "))
	}
	pRef <- match(ref[pok], tAtoms$name)
	nmaxAtoms <- ceiling(mz/tAtoms$mass[pRef])
	vnum <- NULL
	if (class(clist) == "list") {
		vnum <- unlist(clist)
	}
	else {
		vnum <- clist
	}
	dfres <- data.frame(name = ref, mass = tAtoms$mass[pRef],
						valence = tAtoms$valence[pRef])
	if (!is.null(mz))
		dfres$nums <- apply(rbind(unlist(clist), nmaxAtoms),
							2, min)
	return(dfres)
}

formulaToString <- function(vformula,vnames = NULL){
	if(is.null(vnames))	vnames <- names(vformula)
	if(is.factor(vnames)) vnames <- as.character(vnames)
	vnums <- ifelse(vformula>1,vformula,"")
	vnames <- ifelse(as.character(vformula)=="0","",vnames)
	paste(paste(vnames,vnums,sep=""),collapse="")
}

lossesFormulaGeneration <- function(mzrange,atoms=list("C"=15,"H"=50,"O"=20,"N"=6,"P"=2,"S"=1,"Cl"=1),catoms = 1,maxh=1){
	if(!(catoms %in% 0:2)) stop("catoms should be between 0 and 2.")

	tatoms <- makeDfValidAtoms(atoms)

	sig <- diff(mzrange)/2
	mu <- sig+mzrange[1]


	###Finding the formula.
	l_formula <- findFormula(mu,sig,atoms=atoms,mgraph=TRUE,scomp=TRUE)

	if(length(l_formula)==0) return(list())
	###Finding all the heteroatoms
	heteroAtoms <- as.character(tabAtoms()$name[5:11])

	###Now we get the position on the first found formula
	idxHAtoms <- na.omit(match(heteroAtoms,names(l_formula[[1]])))

	to_rm <- which(sapply(l_formula,function(x){sum(x[idxHAtoms]>0)>maxh}))
	if(length(to_rm)!=0){
		l_formula <- l_formula[-to_rm]
	}

	have_heteroatoms <- sapply(l_formula,function(x){any(x[idxHAtoms]>0)})

	mat_f <- do.call(rbind,l_formula)

	masses <- mat_f %*% tatoms$mass

	res <- formulaExtension(masses,mzrange[2],mat_f, have_heteroatoms,catoms )

	oh <- order(res$masses,decreasing=FALSE)
	res$masses <- res$masses[oh]
	res$formula <- res$formula[oh,]
	colnames(res$formula)<-tatoms$name

	return(res)
}
