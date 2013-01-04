##################################################################
##' Map ENTREZ IDs to CNV probes
##'
##' Given a list of ENTREZ IDs, map them to the corresponding CNV
##' probes.
##' @title mapGeneToProbes
##' @param entrezid character vector containing the ENTREZ IDs
##' @return a character vector containing the CNV probes.
##' @author Giovanni d'Ario
##' @export
mapGeneToProbes <- function(entrezid) {
    if(!exists("MarkerRanges"))
        data(MarkerRanges)
    
    GeneRanges <- entrezToGRanges(entrezid)
    probes <- getOverlappingRanges(GeneRanges, MarkerRanges)$probe_id
    return(as.character(probes))
}

