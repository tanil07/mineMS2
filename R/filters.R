#' Title
#'
#' @param m2l The ms2Lib object ot be filtered
#' @param filter_function The filtering function this function take as input a pattern and return TRUE or FALSE.
#' @param args 
#'
#' @return The supplementary ms2Lib arguments.
#' @export
#'
#' @examples
#' print("Exmaples to be put here.")
filterPatterns <-
  function(m2l,filter_function = function(x, threshold, num) {
             return(sum(x@occurences[, "coverage"] > threshold) >= num)
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
