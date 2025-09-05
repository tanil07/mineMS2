# Test discretization errors and warnings {{{1
################################################################################
data(m2l)

##Checking parameters##

#### maxOverlap too high
test_that("maxoverlap_below_0.15", {
    expect_warning(mineMS2::discretizeMzDifferences(m2l, maxOverlap = 0.16))
})

### maxFrags must be below 20
test_that("maxfrags_below_20", {
    expect_error(mineMS2::discretizeMzDifferences(m2l, maxFrags = 21))
})

## count : minimum number of spectra in which the label needs to be found
    test_that("count_conditions", {
    # must be above 2
        expect_warning(mineMS2::discretizeMzDifferences(m2l, count = 1))
    # must be positive 
        expect_error(mineMS2::discretizeMzDifferences(m2l, count = -1))
    # must be inferior to the number of spectra
        expect_error(mineMS2::discretizeMzDifferences(m2l, count = m2l@spectra + 1))
    })

## limMzFormula 
test_that("limmzformula_conditions", {
    # must be a list of 2 floats
        expect_error(mineMS2::discretizeMzDifferences(m2l, limMzFormula = 20)) 
        expect_error(mineMS2::discretizeMzDifferences(m2l, limMzFormula = "14 - 200")) 
    # not too wide interval (limMzFormula[1] <= 10 | limMzFormula[2] > 250)
        expect_warning(mineMS2::discretizeMzDifferences(m2l, limMzFormula = c(1,300)))
})
    

##no heteroatoms
test_that("no_heteroatoms", {
    m2l <- mineMS2::discretizeMzDifferences(m2l, heteroAtoms = FALSE)
    expect_equal(names(m2l@atoms), c("C", "H", "O", "N"))
})

