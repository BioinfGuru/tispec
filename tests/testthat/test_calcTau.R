context("Testing calcTau()")
library(tispec)

test_that('Check output of calcTau', {
    
    expect_identical(calcTau(qnExp), tauExp)
    expect_equal(calcTau(qnExp), tauExp)
    
})