#' @name getDist
#' @title Get tissue specific distribution
#' @description 
#' This function counts the total number of tissue specific protein coding, 
#' lincRNA, and miRNA genes for each tissue.
#' @param x ensembl annotated tau expression fraction file 
#' i.e. output from getMart()
#' @param y threshold tau value, range 0-1
#' @return
#' Returns the number of each type of tissue specific gene that 
#' meet the user defined specificity threshold
#' @examples
#' # Get the number of absolutely specific genes (ASGs, tau = 1) in each tissue
#' head(tauAnno)
#' getDist(tauAnno, 1) 
#' # Get the number of highly specific genes (HSGs, tau = 0.85) in each tissue
#' getDist(tauAnno, 0.85) 
#' @export

getDist <- function(x, y){
    tissues <- names(x[,c(-1, -2, -3)])
    dist <- data.frame()
    for (i in tissues){
        asg <- subset(x, x[,i] >= y)
        t <- length(row.names(asg))
        p <- length(row.names(subset(asg, asg$gene_biotype == "protein_coding")))
        m <- length(row.names(subset(asg, asg$gene_biotype == "lincRNA")))
        l <- length(row.names(subset(asg, asg$gene_biotype == "miRNA")))
        newrow = c(p,m,l,t)
        dist <- rbind(dist,newrow)
    }
    colnames(dist) <- c("protein_coding", "miRNA", "lincRNA", "total")
    row.names(dist) <- tissues
    return(dist)
}