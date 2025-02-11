#' Filter the patterns with intensity coverage
#' 
#' A function to filter the patterns according to their coverage of the dataset. 
#' If less than two spectra in the pattern have a coverage above a threshold of 0.2, then the pattern is discarded.
#' @param m2l The ms2Lib object ot be filtered
#' @param filter_function The filtering function this function take as input a pattern and return TRUE or FALSE.
#' @param args A list of supplementary argument to filter_function
#'
#' @return The supplementary ms2Lib arguments (reducedPatterns argument)
#' @export
#'
#' @examples
#' m2l <- filterPatterns(m2l)
#' print(length(mm2Patterns(m2l)))
#' print(length(mm2ReducedPatterns(m2l)))
filterPatterns <-
  function(m2l,filter_function = function(x, threshold, num) {
             return(sum(x@occurrences[, "coverage"] > threshold) >= num)
           },
           args = list(threshold = 0.2, num = 2)) {
    patt <- mm2Patterns(m2l)
    vval <- sapply(patt, function(x, fun, args) {
      do.call(fun, c(x = x, args))
    }, fun = filter_function, args = args)
    
    ###We remove the non matched patterns.
    m2l@reducedPatterns <- which(vval)
    m2l
}
