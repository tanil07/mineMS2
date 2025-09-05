arrangeSpectra <- function(spectra, step)
{
  for(i in seq(1, length(spectra), by=step))
  {
    
    if(i+(step-1) <= length(spectra))
    {
      gridExtra::grid.arrange(grobs = spectra[i:(i+(step-1))], nrow = step)
    }
    else {
      new_step <- (length(spectra) - i + 1)
      if(new_step == 1)
      {
        gridExtra::grid.arrange(grobs = spectra[length(spectra)], nrow = new_step)
        
      }else
      {
        gridExtra::grid.arrange(grobs = spectra[i:length(spectra)], nrow = new_step)
      }
    }
  }
} 

arrangeSpectra2D <- function(spectra, step_row, step_col)
{
  step_tot <- step_row*step_col
  for(i in seq(1, length(spectra), by=step_tot))
  {
    
    if(i+(step_tot-1) <= length(spectra))
    {
      gridExtra::grid.arrange(grobs = spectra[i:(i+(step_tot-1))],
                              nrow = round(step_tot/2))
    }
    else {
      new_step <- (length(spectra) - i + 1)
      if(new_step == 1)
      {
        gridExtra::grid.arrange(grobs = spectra[length(spectra)], nrow = new_step)
        
      }else
      {
        gridExtra::grid.arrange(grobs = spectra[i:length(spectra)],
                                nrow = round(new_step/2))
      }
    }
  }
}

#' Export in PDF the information about a pattern 
#' 
#' Function to export to PDF the information about a pattern, i.e. list of metabolites containing the pattern, graph of the pattern,
#' list of m/z differences with possible formula and spectra of the metabolites (with colored peaks for vertices in the pattern) 
#' 
#' @param m2l ms2Lib object
#' @param id_p id of the pattern
#' @param infos data frame containing the list of metabolites
#' @param df_mass_diff data frame containing the list of m/z differences
#' @param spectra list of plots for the spectra
#' @param save_dir name of the directory to store the pdf files

printPDF <- function(m2l, id_p, infos, df_mass_diff, spectra, save_dir)
{	
  pdf(file.path(save_dir, paste0("pattern_", id_p, ".pdf")))
  
  grid.text("List of metabolites", draw=TRUE, y = 0.9)
  grid.draw(tableGrob(infos, theme=gridExtra::ttheme_default(base_size = 8)))
  
  plot(m2l,id_p) ## plot the graph of the pattern
  grid.newpage()
  
  
  grid.text("m/z differences in the pattern", draw=TRUE, y = 0.9)
  
  maxrow <- 34
  
  ## case where more than maxrow formulas for one m/z difference
  df_mass_diff_formula <- sapply(df_mass_diff$Formula, FUN = function(x)
  {
    nb <- str_count(x, "\n")
    if(!is.na(nb) && nb > maxrow)
    {
      x_split <- str_split(x, "\n")
      x_split <- sapply(seq_along(x_split[[1]]), function(i)
      {
        print(i)
        if(i%%2 == 0)
        {
          return(paste(x_split[[1]][i], "\n", sep=""))
        }
        else {
          return(paste(x_split[[1]][i], ",", sep=""))
        }
      })
      print(x_split)
      x <- paste(x_split, collapse = "")
      print(x)
    }
    return(x)
  })
  
  df_mass_diff$Formula <- df_mass_diff_formula
  nb_new_lines <- sapply(df_mass_diff_formula, FUN = function(x)
  {
    if(is.na(x)) return(1)
    else return(str_count(x, "\n")+1)
  })
  
  cum_sum <- cumsum(nb_new_lines)
  
  npages <- ceiling(sum(nb_new_lines)/maxrow)
  
  idx <- seq(1, ifelse(length(which(cum_sum >= maxrow)) == 0, nrow(df_mass_diff), which(cum_sum >= maxrow)[[1]]-1))
  
  grid.draw(tableGrob(df_mass_diff[idx,],theme=gridExtra::ttheme_default(base_size = 8), rows = NULL))
  #grid.table(res[idx,], rows=NULL)
  
  if(npages > 1)
  {
    # if several pages needed
    for(i in 2:npages)
    {
      cum_sum_temp <- cumsum(nb_new_lines[(tail(idx,1)+1) :length(nb_new_lines)])
      cum_sum <- c(rep(0,tail(idx,1)), cum_sum_temp)
      grid.newpage()
      if(length(which(cum_sum > maxrow)) == 0)
      {
        idx <- seq(tail(idx,1)+1, length(nb_new_lines))
      }else {
        idx <- seq(tail(idx,1)+1, which(cum_sum > maxrow)[[1]]-1)
      }
      grid.draw(tableGrob(df_mass_diff[idx,], theme=gridExtra::ttheme_default(base_size = 8), rows = NULL))
      #grid.table(res[idx,], rows=NULL)
    }
  }
  
  arrangeSpectra(spectra, step=2)
  dev.off()
}

# plotPatterns ----

#' @rdname plotPatterns
#' @export
setMethod("plotPatterns", "ms2Lib", function(m2l,
                                             ids = NULL,
                                             components = NULL,
                                             occurrences = TRUE,
                                             full = FALSE,
                                             byPage = 9,
                                             titles = NULL,
                                             path_inchi = NULL,
                                             ggplot.l = TRUE,
                                             infos_col = NULL,
                                             save_dir = "none"){
  
  if (!(is.character(save_dir) && length(save_dir) == 1)) {
    stop("'save_dir' argument should be a character of length 1")
  } else if (save_dir != "none") {
    if (!dir.exists(save_dir)) {
      stop(save_dir, " directory not found")
    } else {
      export_pdf <- TRUE
    }
  } else {
    export_pdf <- FALSE
  }
  
  if(is.null(ids)){
    if(length(mm2ReducedPatterns(m2l))==0){
      if(full){
        ids <- vrange(m2l,"P",reduced=FALSE)
      }else{
        stop("No reduced pattern, to plot the full patterns set the 'full' parameter to TRUE")
      }
    }else{
      ids <- vrange(m2l,"P",reduced=TRUE)
    }
  }else{
    tempids <- sapply(ids,parseId,m2l=m2l,simplify = FALSE)
    seq_type <- sapply(tempids,"[[",i="type")
    if(any(seq_type != "patterns")){
      stop("Invalid ids provided ",ids[which(seq_type!="P")]," only patterns ids should be provided")
    }
  }
  
  if(!is.null(titles) & (length(titles)!= length(ids))){
    stop("titles argument should be a character of the same size than ids.")
  }
  
  if(!is.null(components) & (length(components)!= length(ids))){
    stop("components argument should be a list of the same size than ids.")
  }
  
  
  for(pid in seq_along(ids)){
    id <- ids[pid]
    if(!ggplot.l)
    {
      if(!is.null(titles)){
        if(!is.null(components)){
          layout(matrix(c(1,2,2,2),byrow=TRUE,ncol=4))
          plot(0,xaxt="n",yaxt="n",xlab="",ylab="",ylim=c(0,1),xlim=c(0,1),type="n",bty="n")
          text(x = rep(0,length(components[[pid]])),y=seq(0.1,0.9,length=length(components[[pid]])),
               labels = mm2Ids(m2l)[components[[pid]]],cex = 1.0,pos=4)
          plot(m2l,id,title=titles[pid])
          layout(1)
        }else{
          plot(m2l,id,title=titles[pid])
        }
      }else{
        if(!is.null(components)){
          layout(matrix(c(1,2,2,2),byrow=TRUE,ncol=4))
          opar <- par("mar")
          par(mar=c(0.1,0.1,2.1,0.1))
          plot(0,xaxt="n",yaxt="n",xlab="",ylab="",ylim=c(0,1),xlim=c(0,1),type="n",bty="n")
          text(x = rep(0,length(components[[pid]])),y=seq(0.1,0.9,length=length(components[[pid]])),
               labels = mm2Ids(m2l)[components[[pid]]],cex = 1.0,pos=4)
          plot(m2l,id)
          layout(1)
          par(mar=opar)
        }else{
          plot(m2l,id)
        }
      }
      plotOccurrences(m2l,id)
    }
    else
    {	
      p <- support_ggplot(m2l, id, path_inchi = path_inchi, titles = titles)
      infos <- infoPatterns(m2l, id)
      if(!is.null(infos_col))
      {
        infos_to_print <- infos[,infos_col]
      }
      else {
        infos_to_print <- infos
      }			
      res <- listMzDiffbyPattern(m2l, m2l[id], golden_rule=TRUE, spec2Annot=FALSE, export_pdf = export_pdf)
      if(export_pdf)
      {
        printPDF(m2l, id, infos_to_print, res, p, save_dir)
      }
      else
      {
        print(kable(infos_to_print))
        plot(m2l, id)
        print(kable(res))
        #print(p)
        arrangeSpectra2D(p, step_row = 2, step_col = 2)
      }
    }
  }
})
