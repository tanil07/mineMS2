data_dir.c <- system.file("dataset", package = "mineMS2")

NB_SPECTRA <- 51 ## P. nordicum

test_processing_with_metadata <- function() {
    
    path_mgf <- file.path(data_dir.c, 'pnordicum_ms2_spectra.mgf')
    path_metadata <- file.path(data_dir.c, "pnordicum_ms2_info.tsv")
    path_input_graph <- file.path(data_dir.c, 'pnordicum_ms2_gnps.graphml')

    metadata.df <- read.table(path_metadata,
                                    header = TRUE,
                                    sep = "\t") 

    m2l_meta <- mineMS2::ms2Lib(path_mgf, suppInfos = metadata.df)
    testthat::expect_is(m2l_meta, 'ms2Lib')
    
    infos <- mineMS2::getInfo(m2l_meta, "S")
    testthat::expect_is(infos, 'data.frame')
    testthat::expect_true(nrow(infos) == NB_SPECTRA)

    ids <- paste(paste("MZ", infos[,"mz.precursor"]), sep = "_")
    m2l <- mineMS2::setIds(m2l, ids)
    m2l <- mineMS2::discretizeMzDifferences(m2l, dmz = 0.007, ppm = 15,
                                            maxFrags = 15)
    testthat::expect_is(mm2Dags(m2l), 'list')
    testthat::expect_is(mm2Dags(m2l)[[1]], 'igraph')
    testthat::expect_is(mm2EdgesLabels(m2l), 'data.frame')


    m2l <- mineMS2::mineClosedSubgraphs(m2l,
                                        sizeMin = 1,
                                        count = 2)
    testthat::expect_is(mm2Patterns(m2l), 'list')
    testthat::expect_is(mm2Patterns(m2l)[[1]], 'fragPattern')
    
    
    net_gnps <- igraph::read_graph(path_input_graph, "graphml")
    net_gnps <- igraph::simplify(net_gnps,
                                 remove.multiple = FALSE,
                                 edge.attr.comb = "ignore")

    net_gnps <- igraph::as.undirected(net_gnps, mode = "each")

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

testthat::context('Main tests with metadata.')
testthat::test_that('A test processing with metadata works fine.', test_processing_with_metadata())