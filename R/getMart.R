#' @name getMart
#' @title Get ensembl annotations
#' @description
#' Connects to ensembl archive to get id, name and type for every gene
#' @param x species: 'human', 'macaque', or 'mouse'
#' @param y ensembl archive version: 79..90
#' @param z output from calcTau()
#' @return
#' Returns a dataframe containing tau values, tau expression fractions, 
#' ensembl gene ids, names, and types in the selected species from the 
#' selected ensembl archive version.
#' @examples
#' head(tauExp)
#' tauAnno <- getMart(x = 'mouse', y = 79, z = tauExp)
#' head(tauAnno)
#' @export

getMart <- function(x, y, z){
    ## select species
    if (x == 'mouse') {x <- 'mmusculus_gene_ensembl'}
    else if (x == 'human') {x <- 'hsapiens_gene_ensembl'}
    else if (x == 'macaque') {x <- 'mmulatta_gene_ensembl'}
    else {print('species not supported. choose human, macaque, or mouse.')}
    
    ## select host
    if (y == 79) {version.url <- 'http://mar2015.archive.ensembl.org'}
    else if (y == 80) {version.url <- 'http://may2015.archive.ensembl.org'}
    else if (y == 81) {version.url <- 'http://july2015.archive.ensembl.org'}
    else if (y == 82) {version.url <- 'http://sep2015.archive.ensembl.org'}
    else if (y == 83) {version.url <- 'http://dec2015.archive.ensembl.org'}
    else if (y == 84) {version.url <- 'http://mar2016.archive.ensembl.org'}
    else if (y == 85) {version.url <- 'http://july2016.archive.ensembl.org'}
    else if (y == 86) {version.url <- 'http://oct2016.archive.ensembl.org'}
    else if (y == 87) {version.url <- 'http://dec2016.archive.ensembl.org'}
    else if (y == 88) {version.url <- 'http://mar2017.archive.ensembl.org'}
    else if (y == 89) {version.url <- 'http://may2017.archive.ensembl.org'}
    else if (y == 90) {version.url <- 'http://aug2017.archive.ensembl.org'}
    else {print('version not supported. please choose versions 79-90.')}
    
    ## create mart
    mart <- biomaRt::getBM(
        mart = biomaRt::useMart(biomart = "ensembl", dataset = x, host = version.url),
        attributes = c('ensembl_gene_id', 'external_gene_name', 'gene_biotype')
        )
    
    ## set gene ids as row names
    mart$ensembl_gene_id <- as.character(mart$ensembl_gene_id)
    row.names(mart) <- mart$ensembl_gene_id
    mart <- mart[,c(2:ncol(mart))]
    
    ## filter + reorder mart
    mart <- mart[row.names(mart) %in% row.names(z),]
    mart <- mart[ order(row.names(mart)), ]
    cbind(mart, z)
}