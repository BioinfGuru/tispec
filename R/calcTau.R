#' @name calcTau
#' @title Implements tissue specificity algorithm
#' @description 
#' This function first places all quantile normalised values into decile bins. 
#' For each gene, it then uses the bin profile (row of bin values) to calculate 
#' a single tau value between 0 (non-specific) and 1 (absolutely specific). 
#' Finally, for each gene it uses the tau value + quantile normalised counts 
#' to identify a tau expression fraction between 0 (non-specific) and 1 
#' (absolutely specific) which indicates the specificity of a single gene for 
#' a single tissue
#' @param  x quantile normalised data frame i.e. output from quantNorm()
#' @return 
#' Returns a data frame of tau expression fractions indicating the specificty 
#' of a gene for a tissue between 0 (not specific) and 1 (absolutely specific) 
#' @examples
#' head(qnExp)
#' tauExp <- calcTau(qnExp)
#' head(tauExp)
#' @export

calcTau <- function(x){
    xRows <- row.names(x)
    bin <- sapply(x, function(x){
        y <- x[x != 0]
        decile.limits <- stats::quantile(y, probs = seq(.1,.9,.1))
        x[x == 0] <- NA
        x <- lapply(x, function(x){1 + sum(x > decile.limits)})
        x[is.na(x)] <- 0
        return(as.numeric(x))
    })
    bin <- as.data.frame(bin)
    row.names(bin) <-xRows
    tau <- (rowSums(1-(bin/(do.call(pmax, bin)))))/(length(bin)-1)
    tau <- data.frame(tau = tau)
    z<- merge(tau, x, by = "row.names", all = TRUE)
    z$Row.names <- as.character(z$Row.names)
    row.names(z) <- z$Row.names
    z <- z[,c(2:ncol(z))]
    z <- round(z[,1]*(z[,-1]/do.call(pmax, z[,-1])),digits=3)
    z<- merge(round(tau, digits = 3), z, by = "row.names", all = TRUE)
    z$Row.names <- as.character(z$Row.names)
    row.names(z) <- z$Row.names
    z <- z[,c(2:ncol(z))]
    return(z)
}
