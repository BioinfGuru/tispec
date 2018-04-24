context("Testing getTissue()")
library(tispec)

test_that('Check output of getTissue', {
    
    expect_identical(getTissue('tissueA', qnExp, tauAnno), tissueA)
    expect_equal(getTissue('tissueA', qnExp, tauAnno), tissueA)
    
})