#' @name getOptimum
#' @title Extract a set of genes that are both highly specific for and highly 
#' expressed in a tissue of interest.
#' @description 
#' This function ranks genes by their score (as calculated in 
#' getTissue()) and extracts the highest scoring genes that are highly specific
#' to the tissue of interest (>=0.85 tau expression fraction). The number of 
#' genes extracted is defined by the user (z). This function is ideal for 
#' selecting a gene set that has little expression outside the tissue of 
#' interest, while also enough expression in the tissue of interest to 
#' facilitate bench work in the laboratory.
#' @param x output from getMart()
#' @param y output from getTissue()
#' @param z number of genes
#' @return
#' Returns a list of 2 objects: a dataframe and a ggplot barplot. The dataframe is a 
#' subset of the output of getTissue(), containing only the optimum gene set.
#' The barplot shows only the tissues in which the optimum genes are expressed 
#' and the specificity of those genes for those tissues. Red dotted line: 0.85 
#' tau expression fraction threshold for high tissue specificity
#' @examples
#' optimum <- getOptimum(tauAnno, tissueA, 5)
#' optimum$dataframe
#' optimum$barplot
#' @import tidyverse
#' @export
getOptimum <- function(x, y, z){
    
    # ranking
    y$score <- as.numeric(as.character(y$score))
    y$tau <- as.numeric(as.character(y$tau))
    y <- y[order(-y$score, -y$tau), ]
    y <- y[y$frac >= 0.85, ]
    top <- utils::head(y, n = z)
    
    # add the other tissues
    temp <- x[rownames(top), c(-1, -2, -3)]
    temp <- temp[, colSums(temp) != 0]
    temp <- cbind(top, temp)
    temp <- temp[, -c(2:6)]
    temp <- tidyr::gather(temp, 'tissue', 'frac', -'external_gene_name')
    
    # plot
    plot <- ggplot2::ggplot(temp, ggplot2::aes_(x = quote(tissue), y = quote(frac), fill = quote(external_gene_name)))+
        ggplot2::geom_bar(position = "dodge",stat = "identity",colour = "black", width = 0.5)+
        ggplot2::geom_hline(ggplot2::aes(yintercept = 0.85), linetype = "dotted", colour = "red")+
        ggplot2::labs(title = paste('Optimum Gene Set', sep = ""))+
        ggplot2::scale_y_continuous(name = 'Tau Expression Fraction', limits = c(0,1), breaks = round(seq(min(0), max(1), by = 0.1),1))+
        ggplot2::theme_bw()+
        ggplot2::theme(
            legend.title = ggplot2::element_blank(),
            legend.text = ggplot2::element_text(size = 12),
            plot.title = ggplot2::element_text(size = 20, face = "bold", hjust = 0.5),
            axis.title.x = ggplot2::element_blank(),
            axis.text.x = ggplot2::element_text(size = 18, angle = 35, hjust = 1, face = "bold"),
            axis.title.y = ggplot2::element_text(size = 18, face = "bold"),
            axis.text.y = ggplot2::element_text(size = 12),
            panel.border = ggplot2::element_rect(colour = "BLACK",size = 0.5),
            panel.grid.major = ggplot2::element_blank(),
            panel.grid.minor = ggplot2::element_blank(),
            strip.background = ggplot2::element_rect()
        )
    
    ## return a list of objects
    results <- list(top, plot)
    names(results) <- c('dataframe', 'barplot')
    return(results)
}