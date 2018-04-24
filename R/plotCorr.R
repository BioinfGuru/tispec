#' @name plotCorr
#' @title Plot correlation between expression and specificity
#' @description
#' This function plots quantile normalised expression values versus
#' tissue specificity values of all genes in all tissues to visualise 
#' correlation
#' @param x output from getTissue()
#' @return
#' Returns a list object containing 2 ggplot objects: global, local. Global is 
#' the result of plotting QN expression v the overall specificity of each gene 
#' for ANY tissue (tau). Local is the result of plotting QN expression v the 
#' specificity of each gene for the tissue being plotted (tau expression 
#' fraction).
#' 
#' Yellow: tau expression fraction >= 0.85. 
#' 
#' Orange: tau expression fraction = 1. 
#' 
#' Red: generalised additive model trend line
#' 
#' Swatch: 0.99 confidence interval
#'
#' r: correlation coefficient  
#' @examples
#' corrPlots <- plotCorr(tissueA)
#' #corrPlots$global
#' #corrPlots$local
#' @import tidyverse
#' @export

plotCorr <- function(x){
    
    ## QN Exression v Global Specificity (tau)
    asg <- x[x$frac == 1,]
    #hsg <- x[x$frac >= 0.85,]
    cor <- stats::cor.test(x$tau, x$qn, method = c("pearson"))
    corEst <-cor$estimate                  
    corEst <- paste("r = ",round(corEst, digits = 2),sep = "")
    tauPlot <- ggplot2::ggplot(x, ggplot2::aes_(x = quote(tau), y = quote(qn)))+
        ggplot2::geom_point(colour = '#56B4E9' , size = 1, alpha = 0.4)+
        #ggplot2::geom_point(data = hsg, colour = '#F0E442', size = 1.5, alpha = 0.8)+
        ggplot2::geom_point(data = asg, colour = '#E69F00', size = 1.5, alpha = 0.8)+
        ggplot2::guides(fill = FALSE)+
        ggplot2::annotate("text", x = 0.95, y = 14, label = corEst, size = 6)+
        ggplot2::xlab('Specificity (tau)')+
        ggplot2::ylab('QN Expression')+
        ggplot2::scale_x_continuous(limits = c(0,1), breaks = round(seq(min(0), max(1), by = 0.2), 1))+
        ggplot2::scale_y_continuous(limits = c(0,14), breaks = round(seq(min(0), max(14), by = 2), 1))+
        ggplot2::stat_smooth(method = 'auto', level = 0.99, colour = "red", linetype = "dashed", size = 1, na.rm = TRUE)+
        ggplot2::theme_bw()+
        ggplot2::theme(
            plot.title = ggplot2::element_text(size = 20, face = "bold", hjust = 0.5),
            axis.title.x = ggplot2::element_text(size = 25, face = "bold"),
            axis.title.y = ggplot2::element_text(size = 25, face = "bold"),
            axis.text.x = ggplot2::element_text(size = 12, angle = 0, hjust = 0.5), 
            axis.text.y = ggplot2::element_text(size = 12),
            panel.border = ggplot2::element_rect(colour = "BLACK",size = 0.5),
            panel.grid.major = ggplot2::element_blank(),
            panel.grid.minor = ggplot2::element_blank(),
            strip.background = ggplot2::element_rect()
        )
    
    ## QN Expression v Local Specificity (frac)
    asg <- x[x$frac == 1,]
    hsg <- x[x$frac >= 0.85,]
    cor <- stats::cor.test(x$frac, x$qn, method = c("pearson"))
    corEst <-cor$estimate
    corEst <- paste("r = ",round(corEst, digits = 2),sep = "")
    fracPlot <- ggplot2::ggplot(x, ggplot2::aes_(x = quote(frac), y = quote(qn)))+
        ggplot2::geom_point(colour = '#56B4E9' , size = 1, alpha = 0.4)+
        ggplot2::geom_point(data = hsg, colour = '#F0E442', size = 1.5, alpha = 0.8)+
        ggplot2::geom_point(data = asg, colour = '#E69F00', size = 1.5, alpha = 0.8)+
        ggplot2::guides(fill = FALSE)+
        ggplot2::annotate("text", x = 0.95, y = 14, label = corEst, size = 6)+
        ggplot2::xlab('Specificity (frac)')+
        ggplot2::ylab('QN Expression')+
        ggplot2::scale_x_continuous(limits = c(0,1), breaks = round(seq(min(0), max(1), by = 0.2), 1))+
        ggplot2::scale_y_continuous(limits = c(0,14), breaks = round(seq(min(0), max(14), by = 2), 1))+
        ggplot2::stat_smooth(method = 'auto', level = 0.99, colour = "red", linetype = "dashed", size = 1, na.rm = TRUE)+
        ggplot2::theme_bw()+
        ggplot2::theme(
            plot.title = ggplot2::element_text(size = 20, face = "bold", hjust = 0.5),
            axis.title.x = ggplot2::element_text(size = 25, face = "bold"),
            axis.title.y = ggplot2::element_text(size = 25, face = "bold"),
            axis.text.x = ggplot2::element_text(size = 12, angle = 0, hjust = 0.5), 
            axis.text.y = ggplot2::element_text(size = 12),
            panel.border = ggplot2::element_rect(colour = "BLACK",size = 0.5),
            panel.grid.major = ggplot2::element_blank(),
            panel.grid.minor = ggplot2::element_blank(),
            strip.background = ggplot2::element_rect()
        )
    ## return a list of objects
    results <- list(tauPlot, fracPlot)
    names(results) <- c('global', 'local')
    return(results)
}