#' Export Losses Table
#'
#' Export the losses table aith the full informations
#' @param m2l An ms2Lib object
#'
#' @return A data.frame with the losses
#' @export
#'
#' @examples
#' #Loading the data
#' data(m2l)
#' 
#' head(lossesTable(m2l))
lossesTable <- function(m2l){
	if(nrow(mm2EdgesLabels(m2l))==0){
		stop("No loss labels found, use the 'discretizeMassLosses' loss functions")
	}
	return(mm2EdgesLabels(m2l)[,c("lab","mz","mzmin","mzmax","count","formula"),drop=FALSE])

}

#' Export Losses Table with information about the patterns containing those losses
#'
#' 
#' @param m2l An ms2Lib object
#'
#' @return A data.frame with the losses
#' @export
#'
lossesTableComplete <- function(m2l){
	if(nrow(mm2EdgesLabels(m2l))==0){
		stop("No loss labels found, use the 'discretizeMassLosses' loss functions")
	}
    if(length(mm2Patterns(m2l)) == 0)
    {
        stop("No patterns")
    }

    all_pertes <- mm2EdgesLabels(m2l)[,c("lab","mz","mzmin","mzmax","count","formula"),drop=FALSE]
	groups_patterns <- rep("", length(all_pertes[,"lab"]))
    nb_patterns <- rep(0, length(all_pertes[,"lab"]))

    patterns <- mm2Patterns(m2l) ## all patterns

    i = 1
    for(p in patterns)
    {
        print(i)
        g <- mm2Graph(p) ## graph of a pattern
        num_pertes <- unique(edge_attr(g, name="lab"))
        groups_patterns[num_pertes] <- paste(groups_patterns[num_pertes], i, sep=",")
        nb_patterns[num_pertes] <- nb_patterns[num_pertes] + 1
        i <- i+1
    }

    all_pertes[, "nb patterns"] <- nb_patterns
    all_pertes[, "patterns"] <- groups_patterns
    return(all_pertes)
}


formulaSpec2Annot <- function(mass, golden_rule, atoms, index) 
{ 
    #print(mass)
    annot <- find_compo_from_mass(
        mass_target = mass,
        ppm = 15,
        use_golden_ratio = golden_rule,
        elements_vc = NULL)
    #print(annot)
    formula <- apply(annot, 1, function(x){
    res <- sapply(seq_along(x)[1:(length(x)-3)], function(y){
                    if(x[y] == 0) return("")
                    else if (x[y] == 1) return(atoms[y]) 
                    else return(paste(atoms[y], x[y],sep=""))
                })
                res_string <- paste(res[order(index)], collapse="")
                if(grepl("Br|K|Si|Na|F", res_string)) return(NULL)
                return(res_string)
            })
    #print(formula)
    formula <- formula[lapply(formula, length) > 0] ## pour enlever les NULL
    #print(formula)
    return(formula)
}

#' Export Losses Table with information about one pattern
#' mean mz values of losses for spectra of the pattern
#' min, mean and max mz values of losses calculated with all the spectra containing the loss
#' formula calculated with the inner algorithm and the package Spec2Annot
#' 
#' @param m2l An ms2Lib object
#'
#' @return A data.frame with the losses
#' @export
#'
listLossbyPattern <- function(m2l, pattern, golden_rule, spec2Annot, export_pdf)
{
    graph_pattern <- mm2Graph(pattern)
    elabs <- mm2EdgesLabels(m2l)

    atoms <- names(mm2Atoms(m2l))

    
    ## list of mass differences (mean of mz values of the losses in the molecules of the pattern)
    massdiff <- getMzDiff(graph_pattern,mm2Occurences(pattern),mm2Dags(m2l),elabs$mz)
    massdiff <- apply(massdiff,2,mean)

    massdiff_ref <- elabs[edge_attr(graph_pattern, name="lab"),c("mz", "mzmin", "mzmax")]
    massdiff_ref <- apply(massdiff_ref, 2, function(x){
                return(sapply(x, round, digits=4))
    })

    ## if we want to have the chosen formula 
    toccs <- pattern@occurences[,1]
    formula <- get_formula(m2l)[toccs]
    
    oformula <- formula[mm2Occurences(pattern)[,"gid"]]
	filled_formula <- which(sapply(oformula,function(x){(length(x)>0)})&
		  	                          sapply(oformula,function(x){!is.na(x)}))
	if(length(filled_formula)==0){
		  	oformula <- list()
	}else{
		  	oformula <- oformula[filled_formula]
	}
    
    annot <- annotateAFG(pattern,atoms,mm2Dags(m2l),elabs,oformula = oformula)

    chosen_formula <- sapply(annot$edges$formula,function(x,def){
    if(length(x)==0) return(def)
    as.character(x)
  },def=HIGH_MASS_LEGEND)

    ## list of all possible formula ordered by ascending mass
    all_edges_lab <- edge_attr(graph_pattern,name = "lab")
    allf_and_chosen <- elabs$formula[all_edges_lab]
    list_loss_formula <- sapply(elabs$formula[all_edges_lab],function(x,atoms){
      LossFormulaFromSingleString(x,ref = atoms,sep = "|")
    },atoms=atoms,simplify=TRUE)

    m_atoms <- getAtomMass(atoms)
    ## formula sorted by mass
    allf <- sapply(seq_along(allf_and_chosen), function(x){
        loss_formula <- list_loss_formula[[x]]
        formula <- sapply(seq_along(loss_formula@formula[,1]), function(y){formulaToString(loss_formula@formula[y,])})
        if(length(formula) == 0) return(NA)
        mass <- as.numeric(as.character(unlist(loss_formula@formula %*% m_atoms)))
        df <- data.frame(formula=formula, mass=mass)
        df <- df[order(df$mass),]
        if(export_pdf)
        {
            return(paste(df[,'formula'], collapse="\n"))
        }
        else {
            return(paste(df[,'formula'], collapse=","))
        }
        
    }, simplify = TRUE)

    ## spec2Annot
    if(spec2Annot)
    {
        if(golden_rule) {
            atoms_annot <- c('Br', 'Cl', 'S', 'P', 'Si', 'F', 'O', 'N', 'C', 'H')
            index <- c(7,6,5,8,9,10,4,3,1,2)
        }
        else {
        atoms_annot <- c('Br', 'K', 'Cl', 'S', 'P', 'Si', 'Na', 'F', 'O', 'N', 'C', 'H')
        index <- c(7,11,6,5,8,9,12,10,4,3,1,2)
        }

        df <- data.frame()
        for(i in seq_along(massdiff))
        {
            formula <- formulaSpec2Annot(massdiff[i], golden_rule, atoms_annot, index) 
            formula_ref <- formulaSpec2Annot(massdiff_ref[i, "mz"], golden_rule, atoms_annot, index)
            row <- c(massdiff[i], paste(formula, collapse="\n"), paste(formula_ref, collapse="\n"))
            df <- rbind(df, row)
        }
        colnames(df) <- c("mz", "possible_formula_mean_pat", "possible_formula_mean_tot")
        if(export_pdf)
        {
            res <- data.frame("Mass\ndiff\nmean" = massdiff_ref[,"mz"], 
                                "Mass\ndiff\nmin" = massdiff_ref[,"mzmin"], 
                                "Mass\ndiff\nmax" = massdiff_ref[,"mzmax"], 
                                "Formula" = allf, "Chosen\nformula" = chosen_formula, 
                                "Formula\nspec2Annot\nmean" = df[,"possible_formula_mean_tot"])
            res <- res[order(res$"Mass\ndiff\nmean\npat"),]
        }
        else {
            res <- data.frame( "Mass_diff_mean" = massdiff_ref[,"mz"], 
                                "Mass_diff_min" = massdiff_ref[,"mzmin"], 
                                "Mass_diff_max" = massdiff_ref[,"mzmax"], 
                                "Formula" = allf, "Chosen_formula" = chosen_formula, 
                                "Formula_spec2Annot_mean" = df[,"possible_formula_mean_tot"])
            res <- res[order(res$"Mass_diff_mean_pat"),]
        }
    }
    else
    {
        if(export_pdf)
        {
            res <- data.frame("Mass\ndiff\nmean" = massdiff_ref[,"mz"], 
                                "Mass\ndiff\nmin" = massdiff_ref[,"mzmin"], 
                                "Mass\ndiff\nmax" = massdiff_ref[,"mzmax"], 
                                "Formula" = allf, 
                                "Chosen\nformula" = chosen_formula)
            res <- res[order(res$Mass.diff.mean),]
        }
        else {
           res <- data.frame("Mass_diff_mean" = massdiff_ref[,"mz"], 
                                "Mass_diff_min" = massdiff_ref[,"mzmin"], 
                                "Mass_diff_max" = massdiff_ref[,"mzmax"], 
                                "Formula" = allf, 
                                "Chosen_formula" = chosen_formula)
            res <- res[order(res$"Mass_diff_mean"),]
        }
    }
    return(res[!duplicated(res), ])
    #return(list("massdiff" = massdiff, "all_formula" = list(allf), "chosen_formula" = list(chosen_formula)))
}