# Test export infos table {{{1
################################################################################

data(m2l)
data_dir.c <- system.file("dataset", package = "mineMS2")

test_that("export_mzdifftable",{
    expect_no_error(dfDiff <- mzDiffTable(m2l))
    expect_is(dfDiff, 'data.frame')
    expect_true(nrow(dfDiff) == nrow(mm2EdgesLabels(m2l)))

    path_mgf <- file.path(data_dir.c, 'pnordicum_ms2_spectra.mgf')
    m2l_without_edges <- mineMS2::ms2Lib(path_mgf)
    expect_error(mzDiffTable(m2l_without_edges))
})

test_that("export_mzdifftablewithpatterns",{
    expect_no_error(dfDiff <- mzDiffTableComplete(m2l))
    expect_is(dfDiff, 'data.frame')
    expect_true(nrow(dfDiff) == nrow(mm2EdgesLabels(m2l)))

    path_mgf <- file.path(data_dir.c, 'pnordicum_ms2_spectra.mgf')
    m2l_without_edges <- mineMS2::ms2Lib(path_mgf)
    expect_error(mzDiffTableComplete(m2l_without_edges))

    m2l_without_patterns <- mineMS2::discretizeMzDifferences(m2l_without_edges, dmz = 0.007, ppm = 15,
                                            maxFrags = 15)
    expect_error(mzDiffTableComplete(m2l_without_patterns))
})
