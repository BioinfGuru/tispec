#' @name log2Tran
#' @title Log Transforms a data frame
#' @description 
#' Adds a pseudocount (1e-10) to remove zeros, then log2 transforms each 
#' column independently. This normalises the data frame within each column 
#' only (not across columns - see tispec::QuantNorm). After log2 
#' transformation, it resets all negative values to zero, thus setting a 
#' noise threshold of log2(1). Finally it removes all rows where every column 
#' is 0 (i.e. below threshold, removing genes not expressed in any tissue)
#' @param x data frame of mean expression (cells) of genes (rows) in 
#' tissues (columns)
#' @return
#' Returns a data frame of log2 normalised counts, with non-expressed genes 
#' (rows) removed
#' @examples
#' head(meanExp)
#' log2Exp <- log2Tran(meanExp)
#' head(log2Exp)
#' @export

log2Tran <- function(x){
    x[x == 0] <- 0.0000000001
    y <- log2(x)
    y[y<0] <- 0
    y <- y[rowSums(y) != 0,]
    return(y)
}