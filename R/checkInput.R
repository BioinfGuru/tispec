#' @name checkInput
#' @title Checks the input data frame for incompatible values
#' values, and infinite values
#' @description 
#' The input dataframe contains expression data so each cell should contain only 
#' non-negative, non-infinite, numerical characters (>=0). This function 
#' identifies if the dataframe contains character strings, NAs, negative values 
#' or infinite values that must first be removed by the user.
#' @param x data frame of mean expression (cells) of genes (rows) in 
#' tissues (columns)
#' @return Returns errors when incompatible values are found. Also returns a 
#' printout of the number of errors in each column. Nothing is returned If the 
#' dataframe is compatible with the package
#' @examples
#' ## check example data
#' checkInput(meanExp)
#' 
#' ## check data containing character strings
#' #df <- data.frame('col1' = c(1,'foo',3), 'col2' = c(4,5,'bar'))
#' #checkInput(df)
#' 
#' ## check data containing other incompatible values
#' #df <- data.frame('col1' = c(1,2,Inf), 'col2' = c(-4,-Inf,6))
#' #checkInput(df)
#' @export

checkInput <- function(x){
    y <- sapply(x, function(x){
        if(class(x)=='factor'){stop('Some cells contain character strings. Please find and remove them before continuing.', call. = FALSE)} # catches characters
        x[is.infinite(x) == TRUE] <- NA # catches infinites
        suppressWarnings(x[x < 0] <- NA) # catches negatives
        sum(is.na(suppressWarnings(as.numeric(as.character(x))))) # catches non numerics and NAs
    })
    print(y)
    if(sum(y)>0){stop("Some cells contain NAs, negative values, or infinite values. Please find and remove them before continuing.", call. = FALSE)}
}