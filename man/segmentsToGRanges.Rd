\name{segmentsToGRanges}
\alias{segmentsToGRanges}
\title{segmentsToGRanges}
\usage{
  segmentsToGRanges(segments)
}
\arguments{
  \item{segments}{a data frame containing the segments
  obtained from GenomicRanges, or a subset of them.}
}
\value{
  an object of class \code{GenomicRanges} containing the
  genomic coordinates of the segments.
}
\description{
  Transform segments into \code{GenomicRanges} objects
}
\details{
  This function transforms a data frame containing the
  segments obtained from DNAcopy into objects of class
  \code{GenomicRanges}. Note that chromosomes X and Y are
  represented as numbers (23 and 24, respectively).
}
\author{
  Giovanni d'Ario
}

