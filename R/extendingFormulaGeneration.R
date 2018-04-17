RformulaExtension <- function(mzrange,atoms=list("C"=15,"H"=50,"O"=20,"N"=6,"P"=2,"S"=1,"Cl"=1),catoms = 1,maxh=1){
	if(!(catoms %in% 0:2)) stop("catoms should be between 0 and 2.")

	tatoms <- MS2process:::makeDfValidAtoms(atoms)

	sig <- diff(mzrange)/2
	mu <- sig+mzrange[1]


	###Finding the formula.
	l_formula <- MS2process::findFormulaCpp(mu,sig,atoms=atoms,mgraph=TRUE,scomp=TRUE)

	if(length(l_formula)==0) return(list())
	###Finding all the heteroatoms
	heteroAtoms <- as.character(MS2process:::tabAtoms()$name[5:11])

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
