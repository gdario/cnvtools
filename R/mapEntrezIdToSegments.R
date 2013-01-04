##' Map a set of ENTREZ IDs to the segments contained in a 'segment
##' file'.
##'
##' Given a set of ENTREZ IDs, find the segments in a segmented data set
##' overlapping with the corresponding genes.
##' @title mapEntrezIdToSegments
##' @param entrezid character vector containing the ENTREZ IDs
##' @param SegmentRanges object of class \code{GRanges} containing
##' the output of a segmentation algorithm, like DNAcopy (the output
##' slot).
##' @return a data frame containing the ENTREZ IDs and the segments that
##' overlap with the genomic coordinates of the genes.
##' @author Giovanni d'Ario
##' @export
mapEntrezIdToSegments <- function(entrezid=NULL,
                                  SegRanges=NULL) {
    require(biomaRt)
    require(GenomicRanges)
    
    ## Create the GRanges object containing the genes of interest
    Ranges <- entrezToGRanges(entrezid = entrezid)

    ## Compute the overlap
    ovl <- findOverlaps(Ranges, SegRanges)

    eid <- Ranges$entrezgene[queryHits(ovl)]
    segs <- as.data.frame(mcols(SegRanges)[subjectHits(ovl), ],
                          stringsAsFactors = FALSE)
    
    out <- data.frame(entrez_id = eid, segs,
                      stringsAsFactors = FALSE)
    
    return(out)
}
