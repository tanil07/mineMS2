#' Function to export to PDF the information about a pattern, i.e. list of metabolites containing the pattern, graph of the pattern,
#' list of mass differences with possible formula and spectra of the metabolites (with colored peaks for vertices in the pattern) 
#' 
#' @param m2l ms2Lib object
#' @param id_p id of the pattern
#' @param infos data frame containing the list of metabolites
#' @param res data frame containing the list of mass differences
#' @param spectra list of plots for the spectra

printPDF <- function(m2l, id_p, infos, df_mass_diff, spectra)
{
				pdf(file.path("patterns", paste("pattern_", id_p, ".pdf", sep="")))

				grid.text("List of metabolites", draw=TRUE, y = 0.9)
				grid.draw(tableGrob(infos, theme=ttheme_default(base_size = 8)))

				plot(m2l,id_p) ## plot the graph of the pattern
				grid.newpage()
				
				
				grid.text("Mass differences in the pattern", draw=TRUE, y = 0.9)

				maxrow <- 35
				nb_new_lines <- sapply(df_mass_diff$Formula, FUN = str_count, pattern= "\n") + 1
				cum_sum <- cumsum(nb_new_lines)
				npages <- ceiling(sum(nb_new_lines)/maxrow)

				idx <- seq(1, which(cum_sum >= maxrow)[[1]]-1)
				grid.draw(tableGrob(df_mass_diff[idx,],theme=ttheme_default(base_size = 8) ))
				#grid.table(res[idx,], rows=NULL)
				# if several pages needed
				for(i in 2:(npages+1))
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
					grid.draw(tableGrob(df_mass_diff[idx,], theme=ttheme_default(base_size = 8)))
					#grid.table(res[idx,], rows=NULL)
				}
				do.call("grid.arrange", c(spectra, ncol=2))
				dev.off()
}

#' Plot patterns and occurences of an ms2Lib object
#'
#' Plot patterns graph and occurences of an ms2Lib object.
#'
#' @param m2l An ms2Lib object
#' @param ids The ids to be plotted or NULL if all the ids are supposed to be plot
#' @param components An id giving the component of each spectrum to be plotted in first page.

#' @param occurences Shall include the occurences be plotted aswell as the graph.
#' @param full Shall the full patterns be plotted (it can take some times)
#' @param byPage Maximum number of occurences by page.
#' @param titles A vector giving the titles of the MS-MS spectra.
#' @param v_ggplot (default TRUE) if TRUE ggplot version for the plotting of the spectra
#' @param export_pdf (for ggplot version) if TRUE export the spectra to PDF
#' 
#' @export
#' 
#' @return Nothing
#'
#' @examples
#' #Loading the data
#' data(m2l)
#' 
#' #plotting the patterns
#' plotPatterns(m2l,c("P30","P51"))
plotPatterns <- function(m2l,ids=NULL,components = NULL,occurences=TRUE,full=FALSE,byPage=9,titles=NULL,v_ggplot=TRUE, export_pdf=FALSE, inchi=NULL){


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
			stop("Invalid ids furnished ",ids[which(seq_type!="P")]," only patterns ids may be furnished")
		}
	}

	if(!is.null(titles) & (length(titles)!= length(ids))){
		stop("titles argument should be a chracter of the same size than ids.")
	}

	if(!is.null(components) & (length(components)!= length(ids))){
		stop("components argument should be a list of the same size than ids.")
	}

	for(pid in seq_along(ids)){
		id <- ids[pid]
		if(!v_ggplot)
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
				plotOccurences(m2l,id)
		}
		else
		{
			if(!is.null(inchi))
			{
				p <- plot_pattern_ggplot(m2l, id, path_inchi = inchi)
			}
			else {
				p <- plot_pattern_ggplot(m2l, id)
			}
			infos <- infoPatterns(m2l, id)
			res <- listLossbyPattern(m2l, m2l[id], golden_rule=TRUE, spec2Annot=FALSE, export_pdf = export_pdf)
			if(export_pdf)
			{
				if(!dir.exists("patterns"))
				{
					dir.create("patterns")
				}
				printPDF(m2l, id, infos, res, p)
			}
			else
			{
				print(infos)
				plot(m2l, id)
				print(res)
				print(p)
			}
		}
	}
}
