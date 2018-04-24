#' @name plotTop
#' @title Plot a set of highly specific, highly expressed genes in for 
#' tissue of interest
#' @description 
#' This function extracts all genes highly specific to the tissue of 
#' interest (tauExpFrac >=0.85) and ranks them by expression
#' It then selects z number of top expressed genes and plots tau expression 
#' fraction of those genes in all tissues they are expressed in
#' @param x output from getMart()
#' @param y output from getTissue()
#' @param z number of genes
#' @return
#' Returns a barplot showing all tissues in which the set of genes are 
#' expressed and the specificity of those genes for those tissues. 
#' Red dotted line: 0.85 tau expression fraction cut off for high tissue 
#' specificity
#' @examples
#' plotTop(tauAnno, tissueA, 5)
#' @import tidyverse
#' @export

plotTop <- function(x, y, z){
    hsg <- y[y$frac >= 0.85,]
    hsgTop <- hsg[order(-hsg$qn),][1:z,]
    hsgTopFrac <- x[x$external_gene_name %in% hsgTop$external_gene_name[1:z],]
    hsgTopFrac <- hsgTopFrac[,c(-1, -2, -3)]
    hsgTopFrac <- hsgTopFrac[, colSums(hsgTopFrac) != 0]
    hsgTopFrac <- merge(hsgTop, hsgTopFrac, by = "row.names")
    hsgTopFrac <- hsgTopFrac[,c(-1, -3, -4, -5, -6)]
    hsgTopFrac <- tidyr::gather(hsgTopFrac, 'tissue', 'frac', -'external_gene_name')
    ggplot2::ggplot(hsgTopFrac, ggplot2::aes_(x = quote(tissue), y = quote(frac), fill = quote(external_gene_name)))+
        ggplot2::geom_bar(position = "dodge",stat = "identity",colour = "black", width = 0.5)+
        ggplot2::geom_hline(ggplot2::aes(yintercept = 0.85), linetype = "dotted", colour = "red")+
        ggplot2::labs(title = paste('Highest expressed tissue specific genes', sep = ""))+
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
}