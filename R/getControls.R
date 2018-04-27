#' @name getControls
#' @title Get a set of genes that can be used to differentiate between 2 tissues
#' @description
#' This function selects 2 lists of genes, 1 for each tissue. Each list includes 
#' genes that are highly specific and expressed in one tissue with <0.1 
#' quantile normalised expression in the other. This is useful for 
#' differentiating between tissue types. For example, for selecting positive & 
#' negative controls for qPCR to determine if a tissue sample has been 
#' contaminated. 
#' @param x primary tissue, output from getTissue()
#' @param y secondary tissue, output from getTissue()
#' @return
#' Returns a list object containing 1 dataframe of genes for each of the 2 
#' tissues. The data frames are ranked first by score and then by quantile 
#' normalised expression
#' @examples
#' # Usage
#' controls <- getControls(tissueA, tissueB)
#' # To confirm the difference in expression between tissues
#' x <- rownames(head(controls$tissueA, n = 1))
#' qnExp[x, ]
#' tissueA[x, ]
#' tissueB[x, ]
#' @export

getControls <-function(x, y){
    
    # primary
    namex <-deparse(substitute(x))
    x <- x[x$frac >= 0.85, ]
    x <- x[!row.names(x) %in% row.names(subset(y, y$qn >= 0.1)),]
    x$score <- as.numeric(as.character(x$score))
    x$frac <- as.numeric(as.character(x$frac))
    primary <- x[order(-x$score, -x$frac), ]
    
    # secondary
    namey <-deparse(substitute(y))
    y <- y[y$frac >= 0.85, ]
    y <- y[!row.names(y) %in% row.names(subset(x, x$qn >= 0.1)),]
    y$score <- as.numeric(as.character(y$score))
    y$frac <- as.numeric(as.character(y$frac))
    secondary <- y[order(-y$score, -y$frac), ]
    
    ## return a list of objects
    results <- list(primary = primary, secondary = secondary)
    names(results) <- c(namex, namey)
    return(results)
}
