# Test plotting {{{1
################################################################################
data(m2l)

data_dir.c <- system.file("dataset", package = "mineMS2")

## plot spectrum
test_that("plot_spectrum",{
    expect_no_error(plot(m2l, "S1"))
})

## plot DAG
test_that("plot_dag", {
    expect_no_error(plot(m2l,"D1", tkplot=FALSE))
    expect_no_error(plot(m2l, "D1", tkplot=TRUE))
})

## plot pattern graph
test_that("plot_pattern_graph", {
    expect_no_error(plot(m2l, "P1", tkplot = FALSE))
})

## plotting losses is not possible
test_that("plot_losses", {
    expect_error(plot(m2l, "L1"))
})

## plotting several features in the same call to plot is not possible
test_that("plot_several_features", {
    expect_warning(plot(m2l, c("S1", "S2")))}
)

## plot patterns with occurrences spectra and information
test_that("plot_patterns", {
    expect_no_error(plotPatterns(m2l, "P1"))
    ##export in pdf
    expect_error(plotPatterns(m2l, "P1", save_dir = c("save", "dir")))
    expect_no_error(plotPatterns(m2l, "P1", save_dir = data_dir.c))
})
