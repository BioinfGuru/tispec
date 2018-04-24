#' @name getTissue
#' @title Get expression and specificity of each gene in a single tissue 
#' @description 
#' This function takes a user defined tissue and extracts the ensembl 
#' annotations, normalised expression values, tau values, and tau expression
#' fractions of all genes in that tissue.
#' @param x tissue of interest (case sensitive)
#' @param y output from quantNorm()
#' @param z output from getMart()
#' @return 
#' Returns a data frame of ensembl annotations, normalised expression 
#' and specificity values
#' @examples
#' x <- head(colnames(meanExp), n = 1)
#' x <- getTissue(x, qnExp, tauAnno)
#' head(x)
#' y <- rownames(head(x, n = 1))
#' qnExp[y,]
#' tauAnno[y,]
#' @export

getTissue <- function(x, y, z){
    y <- y[order(rownames(y)),]
    z <- z[order(rownames(z)),]
    temp <- z[,c("external_gene_name", "gene_biotype", "tau", x)]
    temp <- cbind(temp, y[x])
    colnames(temp)[c(4,5)] <- c('frac', 'qn')
    temp[,c(1,2,5,3,4)]
}
