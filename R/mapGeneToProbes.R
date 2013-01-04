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
mapGeneToProbes <- function(entrezid,
                            markerRangesFile = paste("..",
                                "compare_ge_and_cnv",
                                "RData",
                                "MarkerRanges.RData",
                                sep = .Platform$file.sep)) {
    if(!exists("MarkerRanges"))
        load(markerRangesFile)
    
    GeneRanges <- entrezToGRanges(entrezid)
    probes <- getOverlappingRanges(GeneRanges, MarkerRanges)$probe_id
    return(as.character(probes))
}

