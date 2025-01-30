# Test plotting {{{1
################################################################################

test_plotting <- function() {
    data(m2l)

    testthat::expect_no_warning(plot(m2l, "S1"))
    testthat::expect_no_warning(plot(m2l, "D1"))
    testthat::expect_no_warning(plot(m2l, "D1", tkplot=TRUE))
    testthat::expect_no_warning(plot(m2l, "P1"))
    testthat::expect_no_warning(plot(m2l, "P1", tkplot=TRUE))

    testthat::expect_no_error(plotPatterns(m2l, "P1"))
    
}

################################################################################
testthat::context('Plotting tests.')
testthat::test_that('All plottings work fine.', test_plotting())