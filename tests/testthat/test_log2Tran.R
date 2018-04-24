context("Testing log2Tran()")
library(tispec)

test_that('Check output of log2Tran', {
    
    expect_identical(log2Tran(meanExp), log2Exp)
    expect_equal(log2Tran(meanExp), log2Exp)
    
})