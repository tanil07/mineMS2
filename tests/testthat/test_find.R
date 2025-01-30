# Test find and select {{{1
################################################################################

test_select_data <- function()
{
    data(m2l)

    ## Find spectra or mass differences by mz
    res1 <- findMz(m2l, 391.1398, dmz = 0.01, type="S")
    testthat::expect_true(all(res1 == c("S32", "S33")))
    
    res2 <- findMz(m2l, 391.1398, type="L")
    testthat::expect_true(length(res2) == 0)

    res3 <- findMz(m2l, 147, "L", dmz = 0.1)
    testthat::expect_true(res3 == "L73")

    ## Find spectra present in patterns
    res1 <- select(m2l, "S2", "P")
    testthat::expect_true(res1 == "P42")

    ## Finding spectra containing mass diff is not implemented
    testthat::expect_error(select(m2l, "S2", "L"))

    ##Finding best covering pattern in terms of intensity for a spectrum
    testthat::expect_error(find(m2l, "S7", "P")) ## needs calculateCoverage before

    m2l <- calculateCoverage(m2l)
    res2 <- find(m2l, "S7", "P")
    testthat::expect_true(res2 == "P47")
}


################################################################################
testthat::context('Tests for selection.')
testthat::test_that('Works fine.', test_select_data())