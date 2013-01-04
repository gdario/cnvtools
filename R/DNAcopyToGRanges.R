##' Convert the segments of a \code{DNAcopy} object into \code{GRanges}
##'
##' This function takes an object of class \code{DNAcopy} and transforms
##' it into an object of class \code{GenomicRanges}. Only
##' the segment-level information is retained, corresponding to the
##' \code{output} component of the \code{DNAcopy} object. The length
##' of the single chromosomes is assumed to be based on HG18, and is
##' automatically set if the \code{chromlen} argument is \code{NULL}.
##' If the segments come from a different genome annotation, they must
##' be provided explicitely as a numeric vector ordered in the same order
##' as the chromosomes in the \code{output} component of \code{x}. This
##' may need some visual inspection, in case \code{random} chromosomes
##' exist.
##' @title DNAcopyToGRanges 
##' @param x an object of class \code{DNAcopy}
##' @param chromlen numeric, the lengths of the chromosomes. See details.
##' @return an object of class \code{GRanges} containing the information
##' contained in the \code{output} component of the \code{DNAcopy}
##' object.
##' @author Giovanni d'Ario
##' @export
DNAcopyToGRanges <- function(x, chromlen=NULL) {
    require(GenomicRanges)

    ## Extract the segments from the DNAcopy object
    segs <- x$output
    
    ## If the chromosomes are represented as numbers, turn them
    ## into characters (for compatibility with most databases)
    if(is.numeric(segs$chrom)) {
        segs$chrom <- as.character(segs$chrom)
        segs$chrom <- ifelse(segs$chrom == "23", "X",
                             ifelse(segs$chrom == "24", "Y",
                                    segs$chrom))
    }
    
    ## Create the GRanges object
    gr <- GRanges(seqnames = segs$chrom,
                  ranges = IRanges(start = segs$loc.start,
                      end = segs$loc.end),
                  ID = segs$ID,
                  num.mark = segs$num.mark,
                  seg.mean = segs$seg.mean)
    if(is.null(chromlen))
        seqlengths(gr) <- chromLengthsHg18$length_bp
    else
        seqlengths(gr) <- chromlen

    return(gr)
}
