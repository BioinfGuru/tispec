#' @name plotDist
#' @title Plot number of tissue specific genes per tissue
#' @description 
#' This function plots the number of absolutely specific genes (ASGs, tau = 1)
#' alongside the number of highly specific genes (HSGs, tau >=0.85) in 
#' each tissue
#' @param x ensembl annotated tau expression fraction file 
#' i.e. output from getMart()
#' @return
#' Returns a barplot of tissue specific gene counts per tissue
#' @examples
#' head(tauAnno)
#' plotDist(tauAnno)
#' @import tidyverse
#' @export

plotDist <- function(x){
    
    ## get asg + hsg
    asg <- getDist(x, 1)[,'total']
    hsg <- getDist(x, 0.85)[,'total']
    
    ## create long data frame
    dist <- data.frame(cbind(asg,hsg))
    tissues <- names(x[,c(-1, -2, -3)])
    dist$tissues <- tissues
    dist <- tidyr::gather(dist, 'spec', 'count', -tissues)
    
    ## plot
    #ggplot2::ggplot(dist, ggplot2::aes(x = stats::reorder(tissues,-count), y = count, fill = as.factor(spec)))+
    ggplot2::ggplot(dist, ggplot2::aes_string(x = 'tissues', y = 'count', fill = 'spec'))+
        ggplot2::geom_bar(position = "dodge",stat = "identity",colour = "black")+
        ggplot2::geom_text(ggplot2::aes_string(label = 'count'), size = 3, vjust = -1, position = ggplot2::position_dodge(width = 1))+
        ggplot2::labs(title = "Distribution") +
        ggplot2::scale_y_continuous(name = "Number of Genes", limits = c(0,2000), breaks = round(seq(min(0), max(2000), by = 200),1))+
        ggplot2::scale_fill_manual(values = c("#0072B2", "#009E73"), name = "")+
        ggplot2::theme_bw()+
        ggplot2::theme(
            legend.text = ggplot2::element_text(size = 12),
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