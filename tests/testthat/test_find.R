# Test selection {{{1
################################################################################


data(m2l)

test_that("find_spectra_by_mz",{
        ## Find spectra by mz
        expect_no_error(res_s <- findMz(m2l, 147.0, type = "S", ppm = 8, dmz = 0.2))
        expect_true(res_s == "S3")
})
    
test_that("find_diff_by_mz",{
    ## Find m/z differences by mz
    expect_no_error(res_l <- findMz(m2l, type = "L", 147.05, dmz = 0.1))
    expect_true(res_l == "L187")
})
    

test_that("find_spectra_in_patterns", {
    ## Find spectra present in patterns
    expect_no_error(res_s3 <- select(m2l, "S3", "P"))
    expect_true(res_s3$S3[1] == "P2")
})

test_that("find_mzdiff_in_patterns", {
    ## Find mzdiff present in patterns
    res_l <- findMz(m2l, type = "L", 147.05, dmz = 0.1)
    expect_no_error(res_p <- select(m2l, res_l, "P"))
    expect_true(res_p$L187[1] == "P70")
})

test_that("find_mzdiff_in_spectra", {
    ## Find mzdiff present in spectra
    res_l <- findMz(m2l, type = "L", 147.05, dmz = 0.1)
    expect_no_error(res_p <- select(m2l, res_l, "S"))
    expect_true(res_p$L187[1] == "S17")
})

test_that("select_with_errors", {
    
    res_l <- findMz(m2l, type = "L", 147.05, dmz = 0.1)
    ## cannot find a mz diff in a dag for now 
    expect_error(res_p <- select(m2l, res_l, "D"))
    ## cannot find a spectrum containing a mz diff for now
    expect_error(res_l <- select(m2l, "S3", "L"))
    ## cannot find a pattern containing a spectrum
    expect_error(res_p <- select(m2l, "P1", "S"))
})
    

test_that("find_best_covering_pattern", {
    # Finding best covering pattern in terms of intensity for a spectrum
    expect_error(find(m2l, "S3", "P")) ## needs calculateCoverage before

    m2l <- calculateCoverage(m2l)
    res_s3 <- find(m2l, "S3", "P")
    expect_true(res_s3 == "P9")
})

test_that("find_all_patterns", {
    ## display all patterns by spectrum
   expect_no_error(all_patterns_by_spectrum <- list_all(m2l,"P", "S", reduced=FALSE))
   ## display all patterns by mz diff
   expect_no_error(all_patterns_by_mzdiff <- list_all(m2l,"P", "L", reduced=FALSE))
})
