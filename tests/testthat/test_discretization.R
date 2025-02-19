# Test discretization errors and warnings {{{1
################################################################################
data(m2l)

##Checking parameters##

#### maxOverlap too high
testthat::expect_warning(mineMS2::discretizeMzDifferences(m2l, maxOverlap = 0.16))
### maxFrags must be below 20
testthat::expect_error(mineMS2::discretizeMzDifferences(m2l, maxFrags = 21))
## count must be above 2
testthat::expect_warning(mineMS2::discretizeMzDifferences(m2l, count = 1))

