#' @name getTissue
#' @title Get expression and specificity of each gene in a single tissue 
#' @description 
#' This function takes a user defined tissue and extracts the ensembl 
#' annotations, normalised expression values, tau values, and tau expression
#' fractions of all genes in that tissue. Each gene's score is the sum of it's 
#' tau expression fraction value and it's 0-1 ranged normalised expression 
#' value. The gene with the greatest expression and highest specificity 
#' recieves the highest score.
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
    if(!x %in% colnames(z)){stop('Check the tissue name you provided', call. = FALSE)}
    df <- z[,c("external_gene_name", "gene_biotype", "tau", x)]
    df <- cbind(df, y[x])
    colnames(df)[c(4,5)] <- c('frac', 'qn')
    df[,c(1,2,5,3,4)]
    
    # get score
    range01 <- function(x){(x-min(x))/(max(x)-min(x))}
    df$qnRanged <- range01(df$qn)
    df$score <- round(df$frac + df$qnRanged, digits = 3)
    df$qnRanged <- NULL
    return(df)
}