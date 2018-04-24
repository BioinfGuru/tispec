#' @name plotDensity
#' @title Plot tissue specificity distribution
#' @description 
#' This function plots the tau value of every gene to see which tau values 
#' occur most often
#' @param x data frame of annotated tissue specificity (tau) values 
#' i.e. output from getMart()
#' @return 
#' Returns a tissue specificity density plot
#' @examples
#' head(tauAnno)
#' plotDensity(tauAnno)
#' @import tidyverse
#' @export

plotDensity <- function(x){
    ggplot2::ggplot(subset(x, select = c('tau')), ggplot2::aes_string(x = 'tau'))+
        ggplot2::geom_density(stat = "density",size = 0.5, fill = "#56B4E9")+
        ggplot2::labs(title = "Tissue Specific Gene Density")+
        ggplot2::scale_x_continuous(limits = c(-0.14,1.14), breaks = round(seq(min(0), max(1), by = 0.1),1), name = "Tau")+    
        ggplot2::scale_y_continuous(limits = c(0,5), breaks = round(seq(min(0), max(6), by = 1),1), name = "Density")+
        ggplot2::theme_bw()+
        ggplot2::theme(
            plot.title = ggplot2::element_text(size = 25,face = "bold", hjust = 0.5),
            axis.title.x = ggplot2::element_text(size = 25, face = "bold"),
            axis.title.y = ggplot2::element_text(size = 25, face = "bold"),
            axis.text.x = ggplot2::element_text(size = 18), 
            axis.text.y = ggplot2::element_text(size = 18),
            panel.border = ggplot2::element_rect(colour = "BLACK",size = 0.5),
            panel.grid.major = ggplot2::element_blank(),
            panel.grid.minor = ggplot2::element_blank(),
            strip.background = ggplot2::element_rect()
        )
}

