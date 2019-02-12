#' @include LossFormula.R
#' @include fragPattern.R

####EXAMPLE DATA
# refa <-c("C","H","O","N")
# ef <- mineMS2:::LossFormula(formula=c("CO2H4","C5H10"),ref=refa)
# sf <- c("H5C2O2","C6H11","C4H18O5")


# vae <- function(x,lf){
#   l_atoms <- colnames(x@formula)
#   lf <- mineMS2:::stringToFormula(lf,vnames=l_atoms)
#   vr <-apply(x@formula,1,is_sub_raw,rf2=lf)
#   vr
# }

#Return the formula witht he most labels..
get_best_label <- function(lf,oformula,mzv,atoms){
  

  
  bloss <- NULL
  bsolo <- 0
  if(is.na(oformula)){
    bloss <- 1:length(lf)
  }else{
    ###Majority vote.
    vm <- sapply(oformula,isSubformula,x=lf,simplify=TRUE)
    bf <- apply(vm,1,sum)
    bloss <- which(bf==max(bf))
  }
  
  ###We keep the one the closest to the furnished mass.
  if(length(bloss)>1){
    m_atoms <- getAtomsMass(atoms)
    vmz <- lf@formula[bloss] %*% m_atoms
    bloss <- bloss[which.min(abs(vmz-mzv))]
    bsolo <- as.numeric(as.numeric(bf[bloss])==nrow(lf@formula))
  }
  
  ###We return the subformula with the maximum size
  c(bloss,bsolo)
  
}


make_labels_raw <- function(g,elabs,atoms,attr_name="lab",formula=character(0),unknow_label = NA_character_){
  
  ####We first extract all the edges labels of the gapths
  re <- edge_attr(g,attr_name)
  ure <- unique(re)
  
  ###We extract all the formula.
  tev <- elabs[ure,,drop=FALSE]
  to_choose <- which(!is.na(tev$formula))

  
  labs <- rep(unknow_label,length(ure))
  all_partial <- rep(NA,length(ure))
  all_solo <- rep(NA,length(ure))
  
  if(length(to_choose) != 0){
  ###For all the formula labels we find the best subformula;
  mlabs <- sapply(to_choose,function(x,elabs,oformula,atoms){
    lf <- LossFormulaFromSingleString(elabs$formula[x],ref=atoms,sep="|")
    bpos <- get_best_label(lf,oformula=oformula,mzv=elabs$mz[x],atoms=atoms)
    ####We generate the formula string.
    c(formulaToString(lf@formula[bpos[1],],vnames=colnames(lf@formula)),ifelse(bpos[2],"TRUE","FALSE"),
      as.character(nrow(lf@formula)))
    },oformula=formula,atoms=atoms,elabs=tev)
  
  labs[to_choose] <- mlabs[1,]
  all_partial[to_choose] <- as.logical(mlabs[2,])
  all_solo[to_choose] <- as.logical(mlabs[3,])
  }
  return(list(id=ure,lab=labs,partial = all_partial,solo=all_solo))
}

