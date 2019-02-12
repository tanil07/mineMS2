#' @include references.R
#' @include allGenerics.R
#' @include allClasses.R

REF_MINEMS2_ATOMS_NAMES <- tabAtoms()$name
REF_MINEMS2_HALOGENS <- tabAtoms()$name[tabAtoms()$halogen]


getAtomsMass <- function(atoms){
  temp <- tabAtoms()
  temp$mass[temp$name %in% atoms]
}

calc_mass_raw <- function(row,atoms_m){
  sum(row*atoms_m)
}


LossFormulaFromSingleString <- function(formula,ref,sep="|"){
  allf <- str_split(formula,fixed(sep),simplify=TRUE)[1,]
  LossFormula(allf,ref=ref)
}


LossFormula <- function(formula=NULL,ref=NULL){
  ####We check the class of the arguments
  lf <- new("LossFormula")
  
  if(is.null(formula)|((length(formula)==1) & (is.na(formula)))){
    if(is.null(ref)){
      lf@formula <- matrix(NA_integer_,nrow=0,ncol=0)
    }else{
      lf@formula <- matrix(NA_integer_,nrow=0,ncol=length(ref))
      colnames(lf@formula) <- ref
    }
  }else if(is.character(formula)){
    if(is.null(ref)){
      if(length(formula)==1){
        lf@formula <- t(sapply(formula,stringToFormula))
      }else{
        temp <- sapply(formula,stringToFormula,simplify=FALSE)
        common_labs <- unique(do.call(c,sapply(temp,names,simplify=FALSE)))
        lf@formula <- t(sapply(formula,stringToFormula,vnames=common_labs))
        
      }
    }else{
      lf@formula <- t(sapply(formula,stringToFormula,vnames=ref))
    }
  }else if(is.matrix(formula)){
    
    lf@formula <- formula
    if(!is.null(ref)){
      colnames(lf@formula) <- ref
    }
  }
  ###We always reorder the formula.
  vf <- match(colnames(lf@formula),REF_MINEMS2_ATOMS_NAMES)
  vf <- order(vf,decreasing = FALSE)
  lf@formula <- lf@formula[,vf,drop=FALSE]
  mode(lf@formula) <- "integer"
  row.names(lf@formula) <- NULL
  return(lf)
}


setMethod("as.character","LossFormula",function(x){
  return(apply(x@formula,1,formulaToString))
})

###Return a Loss formula
combineLossFormula <- function(f1,f2,maxf=3){
            ####We combine all the elements 1 by one.
            if(ncol(f1@formula)!=ncol(f2@formula)) stop("Non compatible formula.")
            if((nrow(f1@formula)==0)|(nrow(f2@formula)==0)) return(LossFormula(ref=ncol(f1@formula)))
  
            maxf1 <- min(maxf,nrow(f1@formula))
            maxf2 <- min(maxf,nrow(f2@formula))
            ###product matrix
            tm <- t(apply(f1@formula[1:maxf1,,drop=FALSE],1,function(x,y){
              t(apply(y,1,function(w,z){
                w+z
              },z=x))
            },y=f2@formula[1:maxf2,,drop=FALSE]))
            # tm <- do.call(rbind,tm)
            tempm <- matrix(c(tm),byrow=FALSE,ncol=ncol(f1@formula))
            LossFormula(tempm,ref=colnames(f1@formula))
}



calcRDBE <- function(lf,vmult=NULL){
  ###We match the input into each cathegory
  vrdbe <- rep(1,nrow(lf@formula))
  
  if(is.null(vmult)){
    vmult <- rep(0,ncol(lf@formula))
    
    pC <- which(colnames(lf@formula)=="C")
    if(length(pC)!=0){
      vmult[pC] <- 1.0
    }
    pH <-  which(colnames(lf@formula)=="H")
    if(length(pH)!=0){
      vmult[pH] <- -0.5
    }
    pN <-  which(colnames(lf@formula)=="N")
    if(length(pN)!=0){
      vmult[pN] <- -0.5
    }
    PHalo <- which(colnames(lf@formula) %in% REF_MINEMS2_HALOGENS)
    if(length(PHalo)!=0){
      vmult[PHalo] <- -0.5
    }
  }
  calcRDBE_raw(lf@formula,vmult)
}

calcRDBE_raw<- function(lfm,vmult){
  return((lfm %*% vmult)+1)
}


convertFormula <- function(x,y){
  vm <- match(names(x),colnames(y))
  val <- rep(0,ncol(y))
  val[vm] <- x
  as.integer(val)
}

orderByRDBE <- function(lf){
  rdbe <- calcRDBE(lf)
  lf@formula[order(rdbe,decreasing = TRUE)]
  lf
}


setMethod("%in%",signature=list(x="character",table="LossFormula"),function(x,table){
  vf <- stringToFormula(x,vnames=colnames(table@formula))
  any(apply(table@formula, 1, function(x, want) isTRUE(all.equal(x, want)), vf))
})


setMethod("addFormula",signature = list(x = "LossFormula",lf = "ANY"),function(x,lf){
  l_atoms <- colnames(x@formula)
  
  if(typeof(lf)=="character"){
    lf <- stringToFormula(lf,vnames=l_atoms)
  }else if(typeof(lf)=="numeric"){
    lf <- convertFormula(lf,x@formula)
  }else if(typeof(lf)=="integer"){
    lf <- convertFormula(lf,x@formula)
  }
  x@formula <- cbind(x@formula,lf)
  x
})

is_sub_raw <- function(rf1,rf2){
  all(rf1<=rf2)
}


###Vectorized implementation.
setMethod("isSubformula",signature = list(x="LossFormula",lf = "character"),function(x,lf){
  l_atoms <- colnames(x@formula)
  lf <- stringToFormula(lf,vnames=l_atoms)
  vr <-apply(x@formula,1,is_sub_raw,rf2=lf)
  vr
})


setMethod("isSubformula",signature = list(x="LossFormula",lf="LossFormula"),function(x,lf){
  
  l_atoms <- colnames(object@formula)
  
  tm <- t(apply(x@formula,1,function(x,y,agg){
    res <- apply(y,1,function(w,z){
      is_sub_raw(w,z)
    },z=x)
    
    any(res)
  },y=lf@formula))
  
  any(tm)
})



setMethod("addFormula",signature = list(x = "LossFormula",lf = "LossFormula"),function(x,lf){
  if(any(colnames(lf@formula)!=colnames(x@formula))) stop("Incompatible atoms in formula.")
  l_atoms <- colnames(lf@formula)
  
  idobj <- as.character(lf)
  idx <- as.character(x)
  
  pok <- !(idx %in% idobj)
  
  if(sum(pok)!=0){
    lf@formula <- rbind(lf@formula,x@formula[pok,,drop=FALSE])
  }
  lf
})


setMethod("show","LossFormula",function(object){
  cat("A LossFormula object containing",nrow(object@formula),"formula with atoms",paste(colnames(object@formula),collapse = ","))
})


setMethod(length,"LossFormula",function(x){
  nrow(x@formula)
})
