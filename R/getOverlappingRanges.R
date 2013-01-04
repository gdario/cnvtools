##' Computes the overlap between a pair of \code{GenomicRanges} objects
##'
##' Compute the overlap between a pair of \code{GenomicRanges} objects
##' @title getOverlappingRanges
##' @param query an object of class \code{GenomicRanges}. The 'query'
##' component in the \code{findOverlaps} function call.
##' @param subject an object of class \code{GenomicRanges}. The
##' 'subject'
##' component in the \code{findOverlaps} function call.
##' @return a data frame containing the coordinates of the 'query' and
##' of the
##' overlapping 'subject' ranges.
##' @author Giovanni d'Ario
##' @export
getOverlappingRanges <- function(query, subject) {
    ## Compute the overlap
    ovl <- findOverlaps(query, subject)
    
    query <- query[queryHits(ovl), ]
    subject <- subject[subjectHits(ovl), ]

    out <- cbind(as(query, "data.frame"),
                 as(subject, "data.frame"))

    return(out)
}
