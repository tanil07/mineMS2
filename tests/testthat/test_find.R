test_select_data <- function()
{
    data(m2l)

    ## Find spectra or m/z differences by mz
    res1 <- findMz(m2l, 147.0, type = "S", ppm = 8, dmz = 0.2)
    testthat::expect_true(res1 == "S3")
    
    res2 <- findMz(m2l, type = "L", 147.05, dmz = 0.1)
    testthat::expect_true(res2 == "L187")

    ## Find spectra present in patterns
    res3 <- select(m2l, "S3", "P")
    testthat::expect_true(res3$S3[1] == "P2")

    ## Finding best covering pattern in terms of intensity for a spectrum
    testthat::expect_error(find(m2l, "S3", "P")) ## needs calculateCoverage before

    m2l <- calculateCoverage(m2l)
    res4 <- find(m2l, "S3", "P")
    testthat::expect_true(res4 == "P9")
}

testthat::context('Tests for selection.')
testthat::test_that('Works fine.', test_select_data())