#' @name quantNorm
#' @title Quantile normalise a data frame
#' @description 
#' This function quantile normalises the entire data frame of log normalised 
#' counts, allowing comparisons across tissues. First, the row and column 
#' names are collected and all zeros are changed to NA. After converting to a 
#' matrix, the data frame is quantile normalised. Finally, all NAs are 
#' reverted back to zero, the matrix is converted back to a data frame and 
#' row and column names are re-attached
#' @param x log transformed data frame of mean expression (cells) of 
#' genes (rows) in tissues (columns)
#' @return
#' Returns a data frame of equal dimensions with quantile normalised 
#' counts across all columns
#' @examples
#' head(log2Exp)
#' qnExp <- quantNorm(log2Exp)
#' head(qnExp)
#' @export

quantNorm <- function(x){
    x.cols <- names(x)
    x.rows <- row.names(x)
    x[x == 0] <- NA
    x_m <- as.matrix(x)
    x <- round(preprocessCore::normalize.quantiles(x_m), digits = 3)
    x[is.na(x)] <- 0
    x <-data.frame(x)
    names(x)[c(seq_along(x.cols))] <- x.cols
    row.names(x)[c(seq_along(x.rows))] <- x.rows
    return(x)
}