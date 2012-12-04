##' Transform segments into \code{GenomicRanges} objects
##'
##' This function transforms a data frame containing the segments
##' obtained from DNAcopy into objects of class \code{GenomicRanges}.
##' Note that chromosomes X and Y are represented as numbers (23 and
##' 24, respectively).
##' @title segmentsToGRanges
##' @param segments a data frame containing the segments obtained
##' from GenomicRanges, or a subset of them.
##' @return an object of class \code{GenomicRanges} containing the
##' genomic coordinates of the segments.
##' @author Giovanni d'Ario
segmentsToGRanges <- function(segments) {
    ## Fix the chromosomes (X is represented as 23 and Y as 24,
    ## and not as letters.
    Segments$chrom <- ifelse(Segments$chrom < 23,
                             Segments$chrom,
                             ifelse(Segments$chrom == 23,
                                    "X", "Y"))
    ## Create the GRanges object for the segments
    SegRanges <- with(Segments,
                      GRanges(seqnames = chrom,
                              ranges = IRanges(start = loc.start,
                                  end = loc.end),
                              seg.mean = seg.mean,
                              ID = ID))
    return(SegRanges)
}
