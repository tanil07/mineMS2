###Function to scroe a pattern as extracted from R


sigmoid <- function(x,weigth=3.5){
	return(1/(1+exp(-weigth*(x-0.5))))
}


#' Determine edges scores
#'
#' @param m2l An ms2Lib object
#' @param knownLoss THe pnality to use for know losses.
#' @param selectedLosses The majorition to use for known loss.
#' @param monoAtom THe panality associated to known loss.
#' @param weigth The weight associated to the sigmoid function.
#'
#' @return
#' @export
#'
#' @examples
#' print("exmaples to be put here")
setMethod("scoresLosses","ms2Lib",function(m2l,knownLoss=0.2,selectedLosses=3,monoAtom=-1,weigth=3.5){
	elabs <- mm2EdgesLabels(m2l)

	m_known <- ifelse(elabs$adv_loss,knownLoss,1)

	a_mono <- apply(elabs[,c("carb_only","nitrogen_only","ch_only")],1,
					function(x){
						any(x)
					})
	a_mono <-ifelse(a_mono,monoAtom,0)

	a_mzv <- log(elabs$mz,100)

	tcount <- elabs$count
	tsum <- sum(tcount)

	probs <- sapply(tcount,function(x,tmax){
		(tmax-x)/tmax
	},tmax=tsum)

	minp <- min(probs)
	maxp <- max(probs)

	###Maw is between 0 and 1.
	probs <- (probs-minp)/(maxp-minp)

	scores <-sigmoid(probs,weigth=weigth)*m_known+a_mono+a_mzv

	elabs$score <- scores
	mm2EdgesLabels(m2l) <- elabs
	m2l
})


#' Score the patterns contained in an ms2Lib object
#'
#' Scores the patterns contained in a ms2Lib object using the scoring  of the edges.
#' shloub b e called after
#'
#' @param m2l An ms2Lib object ot be scored.
#'
#' @return the completed ms2Lib object
#' @export
#'
#' @examples
#' print("examples to be put here")
setMethod("scorePatterns","ms2Lib",function(m2l){
	if(!("score" %in% colnames(mm2EdgesLabels(m2l)))){
		stop("'score' of mass losses needs to be calculated before the patterns are scored, used the scoreLosses")
	}
	lat <- mm2Lattice(m2l)
	first_pattern <- max(mm2LatticeIdxItems(m2l))

	vec_scores <- rep(NA_real_,vcount(lat))

	###We get the scores of the label.
	scores_losses <- mm2EdgesLabels(m2l)$score

	###The score is calculated for each pattern.
	vec_scores[seq(first_pattern+1,length(vec_scores))] <- sapply(mm2Patterns(m2l),function(x,vloss){
		nv <- as.integer(vcount(mm2Graph(x)))
		tx <- as_data_frame(mm2Graph(x),what="edges")
		tx$from <- as.integer(tx$from)
		tx$to <- as.integer(tx$to)
		# cat(class(tx),class(vloss),class(nv))
		scorePattern(tx,vloss,nv)
	},vloss=scores_losses)

	lat <- set_vertex_attr(lat,"score",value = vec_scores)
	# print(vertex_attr(lat,"score"))

	mm2Lattice(m2l) <- lat
	m2l
})


#' Reduce the lattice to k nodes.
#'
#' Reduce the lattice to the top-k patterns.
#'
#' @param m2l ms2Lib object.
#' @param k The number of remaining nodes.
#'
#' @return the filled ms2Lib object.
#' @export
#'
#' @examples
#' print("Examples here")
kLatticeReduction <- function(m2l,k){
	if(length(mm2Patterns(m2l))==0){
		stop("No patterns detected or impossible ot reduce the lattice.")
	}
	###The score vector is captured.
	lat <- mm2Lattice(m2l)
	sc <- vertex_attr(lat,"score")
	num_objects <- length(mm2LatticeIdxItems(m2l))

	if((vcount(lat)-num_objects-k)<=0) stop("Impossible to reduce the lattice to ",k,
											" because there is ony ",vcount(lat)-num_objects," patterns.")

	message("Initial calculations")
    n <- vecInitialisation(mm2Patterns(m2l),sc,num_objects)

    df_vertices <- as_data_frame(lat,"vertices")
    df_vertices$name <- as.integer(df_vertices$name)

    df_edges <- as_data_frame(lat,"edges")
    df_edges$from <- as.integer(df_edges$from)
    df_edges$to <- as.integer(df_edges$to)

    sc <- sc[(num_objects+1):length(sc)]
    # browser()
    message("Lattice reductions")
	removed <-  reduceLatticeK_greedy(df_vertices,df_edges, sc,n,k)

	###We get the indices of the conserved patterns
	mm2ReducedPatterns(m2l) <- as.integer(complementIdx(length(mm2Patterns(m2l)),removed)+1)

	m2l
}
