data_dir.c <- system.file("dataset", package = "mineMS2")

NB_SPECTRA <- 51 ## P. nordicum

test_processing <- function() {
    
    path_mgf <- file.path(data_dir.c, 'dda_msms_pnordicum.mgf')
    path_input_graph <- file.path(data_dir.c, 'graph_gnps_pnordicum.graphml')

    testthat::expect_no_error(checkFormat(path_mgf))
    testthat::expect_error(checkFormat(path_input_graph))

    m2l <- mineMS2::ms2Lib(path_mgf)
    testthat::expect_is(m2l, 'ms2Lib')
    
    infos <- mineMS2::getInfo(m2l, "S")
    testthat::expect_is(infos, 'data.frame')
    testthat::expect_true(nrow(infos) == NB_SPECTRA)
    ## without supp data, information about the spectra must be the following three:
    testthat::expect_true(all(c("mz.precursor", "file", "formula") %in% names(infos)))


    testthat::expect_error(m2l <- mineMS2::setIds(m2l, "mz"))
    ids_error <-  paste(paste("S", infos[,"mz.precursor"], sep = "_"))
    testthat::expect_error(m2l <- mineMS2::setIds(m2l, ids_error))
    ids <- paste(paste("MZ", infos[,"mz.precursor"]), sep = "_")
    m2l <- mineMS2::setIds(m2l, ids)
    testthat::expect_is(mm2Ids(m2l), 'character')
    testthat::expect_equal(mm2Ids(m2l), ids)

    m2l <- mineMS2::discretizeMzDifferences(m2l, dmz = 0.007, ppm = 15,
                                            maxFrags = 15)
    testthat::expect_is(mm2Dags(m2l), 'list')
    testthat::expect_is(mm2Dags(m2l)[[1]], 'igraph')
    testthat::expect_is(mm2EdgesLabels(m2l), 'data.frame')

    testthat::expect_is(mm2NodesLabels(m2l), 'data.frame')
    testthat::expect_true(nrow(mm2NodesLabels(m2l)) == 0 && ncol(mm2NodesLabels(m2l)) == 0)

    m2l <- mineMS2::mineClosedSubgraphs(m2l,
                                        sizeMin = 1,
                                        count = 2)
    testthat::expect_is(mm2Patterns(m2l), 'list')
    testthat::expect_is(mm2Patterns(m2l)[[1]], 'fragPattern')

    testthat::expect_no_error(show(m2l))
    
    net_gnps <- igraph::read_graph(path_input_graph, "graphml")
    net_gnps <- igraph::simplify(net_gnps,
                                 remove.multiple = FALSE,
                                 edge.attr.comb = "ignore")

    net_gnps <- igraph::as_undirected(net_gnps, mode = "each")

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
    
}

testthat::context('Main tests.')
testthat::test_that('A test processing works fine.', test_processing())
