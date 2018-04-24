context("Testing getMart()")
library(tispec)

test_that('Check output of getMart', {
    
    expect_identical(getMart(x = 'mouse', y = 79, z = tauExp), tauAnno)
    expect_equal(getMart(x = 'mouse', y = 79, z = tauExp), tauAnno)
    
})