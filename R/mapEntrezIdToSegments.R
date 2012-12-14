##' Map a set of ENTREZ IDs to the segments contained in a 'segment
##' file'.
##'
##' Given a set of ENTREZ IDs, find the segments in a segmented data set
##' overlapping with the corresponding genes.
##' @title mapEntrezIdToSegments
##' @param entrezid character vector containing the ENTREZ IDs
##' @param segmentFile string containing the full path to a file
##' containing
##' the output of a segmentation algorithm, like DNAcopy (the output
##' slot).
##' @return a data frame containing the ENTREZ IDs and the segments that
##' overlap with the genomic coordinates of the genes.
##' @author Giovanni d'Ario
##' @export
mapEntrezIdToSegments <- function(entrezid, segmentFile=NULL) {

    require(biomaRt)
    require(GenomicRanges)

    if(!exists("Segments"))
        if(is.null(segmentFile))
            data(Segments)
        else
            load(file = segmentFile)

    ## Create the GRanges object containing the genes of interest
    Ranges <- entrezToGRanges(entrezid = entrezid)
    SegRanges <- segmentsToGRanges(segments = Segments)    

    ## Compute the overlap
    ovl <- findOverlaps(Ranges, SegRanges)

    eid <- Ranges$entrezgene[queryHits(ovl)]
    segs <- mcols(SegRanges)[subjectHits(ovl), ]

    ## Convert segs into a data frame
    segs <- as(segs, "data.frame")
    
    out <- data.frame(entrez_id = eid, segs,
                      stringsAsFactors = FALSE)
    
    return(out)
}
