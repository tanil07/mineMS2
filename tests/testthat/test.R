# vi: fdm=marker ts=4 et cc=80 tw=80

# Constants {{{1
################################################################

TEST.DIR <- file.path(getwd(), '..')
IN.DIR  <- file.path(TEST.DIR, 'inputs')
OUT.DIR  <- file.path(TEST.DIR, 'outputs')

NB_SPECTRA <- 52 ## P. verrucosum

# Test processing {{{1
################################################################################

test_processing <- function() {
    
    path_mgf <- file.path(IN.DIR, 'test.mgf')
    path_input_graph <- file.path(IN.DIR, 'test.graphml')
    path_expected_output_graph <- file.path(IN.DIR, 'test_annotated.graphml')
    path_output_graph <- file.path(OUT.DIR, 'test_processing.graphml')

    m2l <- mineMS2::ms2Lib(path_mgf)
    testthat::expect_is(m2l, 'ms2Lib')
    
    infos <- mineMS2::getInfo(m2l, "S")
    testthat::expect_is(infos, 'data.frame')
    testthat::expect_true(nrow(infos) == NB_SPECTRA)
    ## without supp data, information about the spectra must be the following three:
    testthat::expect_true(all(c("mz.precursor", "file", "title") %in% names(infos)))

    ids <- paste(paste("MZ", infos[,"mz.precursor"]), sep = "_")
    m2l <- mineMS2::setIds(m2l, ids)
    m2l <- mineMS2::discretizeMzDifferences(m2l,
                                            dmz = 0.007,
                                            ppm = 15,
                                            heteroAtoms = TRUE,
                                            maxFrags = 15)
    testthat::expect_is(mm2Dags(m2l), 'list')
    testthat::expect_is(mm2Dags(m2l)[[1]], 'igraph')
    testthat::expect_is(mm2EdgesLabels(m2l), 'data.frame')


    m2l <- mineMS2::mineClosedSubgraphs(m2l,
                                        sizeMin = 1,
                                        count = 2)
    testthat::expect_is(mm2Patterns(m2l), 'list')
    testthat::expect_is(mm2Patterns(m2l)[[1]], 'mineMS2::fragPattern')
    
    mineMS2::plotPatterns(m2l, full = TRUE)
    
    net_gnps <- igraph::read_graph(path_input_graph, "graphml")
    net_gnps <- igraph::simplify(net_gnps,
                                 remove.multiple = FALSE,
                                 edge.attr.comb = "ignore")
    testthat::expect_is(net_gnps, 'igraph')
    components <- mineMS2::findGNPSComponents(net_gnps,minSize=3,pairThreshold = 0.9)
    
    patterns <- mineMS2::findPatternsExplainingComponents(m2l,
                                                          components,
                                                          metric=c("recall",
                                                                   "precision",
                                                                   "size"))
    testthat::expect_is(patterns, 'list')
    
    annotated_net <- mineMS2::annotateNetwork(components,
                                              net_gnps,
                                              patterns)
    testthat::expect_is(annotated_net, 'igraph')
    
    igraph::write_graph(graph = annotated_net,
                        format = "graphml",
                        file = path_output_graph)
    
}

# Test MGF filename {{{1
################################################################################

test_mgf_filename <- function() {
    
    path_mgf <- file.path(IN.DIR, 'test.mgf')
    mgf_with_no_ext_mgf <- file.path(OUT.DIR, 'test.dat')
    file.copy(path_mgf, mgf_with_no_ext_mgf)
    
    m2l <- mineMS2::ms2Lib(mgf_with_no_ext_mgf)
    testthat::expect_is(m2l, 'ms2Lib')
}

# Main {{{1
################################################################################

testthat::context('Main tests.')
testthat::test_that('A test processing works fine.', test_processing())
testthat::test_that('We can provide any filename for a single MGF file', test_mgf_filename())
