context("Testing quantNorm()")
library(tispec)

test_that('Check output of quantNorm', {
    
    expect_identical(quantNorm(log2Exp) , qnExp)
    expect_equal(quantNorm(log2Exp) , qnExp)
    
})