# Test filtering patterns according to coverage {{{1
################################################################################
data(m2l)

test_that("test_filter", {
    expect_error(m2l <- filterPatterns(m2l))
    expect_true(length(mm2Patterns(m2l)) == length(mm2ReducedPatterns(m2l)))
    m2l <- calculateCoverage(m2l)
    expect_no_error(m2l <- filterPatterns(m2l))
    expect_true(length(mm2ReducedPatterns(m2l)) <= length(mm2Patterns(m2l)))
})