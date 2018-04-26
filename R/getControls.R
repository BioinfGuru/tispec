#' @name getControls
#' @title Get a set of genes that can be used to differentiate between 2 tissues
#' @description
#' This function selects 2 set of genes, 1 for each tissue. Each set includes 
#' genes that are highly specific and expressed in one tissue with <0.1 
#' quantile normalised expression in the other. This is useful for 
#' differentiating between tissue types. For example, for selecting positive & 
#' negative controls for qPCR to determine if a tissue sample has been 
#' contaminated.
#' @param x primary tissue, output from getTissue()
#' @param y secondary tissue, output from getTissue()
#' @return
#' Returns a list object containing 1 dataframe of genes for each of the 2 
#' tissues.
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
    
    ## get primary highly specific genes with low expression in secondary
    primary <- x[x$frac >= 0.85,]
    primary <- primary[!row.names(primary) %in% row.names(subset(y, y$qn >= 0.1)),]
    primary <- primary[order(-primary$qn),]
    
    ## get secondary highly specific genes with low expression in primary
    secondary <- y[y$frac >= 0.85,]
    secondary <- secondary[!row.names(secondary) %in% row.names(subset(x, x$qn >= 0.1)),]
    secondary <- secondary[order(-secondary$qn),]
    
    ## return a list of objects
    results <- list(primary = primary, secondary = secondary)
    namex <-deparse(substitute(x))
    namey <-deparse(substitute(y))
    names(results) <- c(namex, namey)
    return(results)
}
