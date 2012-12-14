##' Create a reverse mapped matrix starting from DNAcopy segments
##'
##' This function simply runs \code{mapSegmentsToMarkers} to all the
##' chromosomes listed in the \code{chrom} argument, and binds all
##' the chromosome-specific results into a matrix.
##' @title createRevMappedCN
##' @param x an object of class \code{DNAcopy}
##' @param chrom a numeric vector containing the chromosomes of interest
##' @return a matrix containing as many columns as are the patients in
##' the dataset, and as many rows as are the probes. Each cell contains,
##' for each patient, the mean segment value for the probe falling inside
##' a segment.
##' @author Giovanni d'Ario
##' @export
createRevMappedCN <- function(x,
                              chrom = 1:22) {
    require(DNAcopy)
    Lst <- list()
    for(chr in chrom)
        Lst[[chr]] <- mapSegmentsToMarkers(x, chrom = chrom)
    out <- do.call("rbind", Lst)
    return(out)
}
