#' @name plotGene
#' @title Plot specificity of a user defined gene for each tissue
#' @description 
#' This function takes a user defined gene (case sensitive) and plots the
#' tau expression fraction of that gene in each tissue. This allows the user 
#' to visualize the tissue/tissues in which the genes has relatively 
#' enriched/depleted expression
#' @param x ensembl annotated tau expression fraction file 
#' i.e. output from getMart()
#' @param y user define gene, case dependent, to see gene available names 
#' view output of getMart()
#' @return
#' Returns a bar plot of tau espression fraction in each tissue
#' @examples
#' x <- head(tauAnno$external_gene_name, n = 1)
#' plotGene(tauAnno, x)
#' y <- rownames(subset(tauAnno, tauAnno$external_gene_name == x))
#' meanExp[y,]
#' qnExp[y,]
#' tauAnno[y,]
#' @import tidyverse
#' @export

plotGene <- function(x, y){
    geneFrac <- data.frame(t(subset(x[,c(-1, -2, -3)], x$external_gene_name == y)))
    geneFrac[,1] <- round(geneFrac[,1], digits = 2)
    geneID <- colnames(geneFrac)
    geneTau <- x[colnames(geneFrac), "tau"]
    if (!y %in% x$external_gene_name){
        stop(paste(y, ' : ', geneID, ' : Gene not found, check expression', sep = ''))
    }
    colnames(geneFrac)[1] <- 'tauExpFrac'
    ggplot2::ggplot(geneFrac, ggplot2::aes_(x = rownames(geneFrac), y = quote(tauExpFrac), fill = quote(tauExpFrac)))+
        ggplot2::geom_bar(position = "dodge",stat = "identity",colour = "black")+
        ggplot2::geom_text(ggplot2::aes_(label = quote(tauExpFrac)), size = 3, vjust = -1, position = ggplot2::position_dodge(width = 1))+
        ggplot2::labs(title = y)+
        ggplot2::scale_y_continuous(name = "Tau Expression Fraction", limits = c(0,1), breaks = round(seq(min(0), max(1), by = 0.1),1))+
        ggplot2::theme_bw()+
        ggplot2::theme(
            legend.text = ggplot2::element_text(size = 12),
            legend.position = "none",
            plot.title = ggplot2::element_text(size = 20, face = "bold", hjust = 0.5),
            axis.title.x = ggplot2::element_blank(),
            axis.text.x = ggplot2::element_text(size = 12, angle = 35, hjust = 1, face = "bold"),
            axis.title.y = ggplot2::element_text(size = 12, face = "bold"),
            axis.text.y = ggplot2::element_text(size = 12),
            panel.border = ggplot2::element_rect(colour = "BLACK",size = 0.5),
            panel.grid.major = ggplot2::element_blank(),
            panel.grid.minor = ggplot2::element_blank(),
            strip.background = ggplot2::element_rect()
        )
}