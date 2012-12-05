\name{getOverlappingRanges}
\alias{getOverlappingRanges}
\title{getOverlappingRanges}
\usage{
  getOverlappingRanges(query, subject)
}
\arguments{
  \item{query}{an object of class \code{GenomicRanges}. The
  'query' component in the \code{findOverlaps} function
  call.}

  \item{subject}{an object of class \code{GenomicRanges}.
  The 'subject' component in the \code{findOverlaps}
  function call.}
}
\value{
  a data frame containing the coordinates of the 'query'
  and of the overlapping 'subject' ranges.
}
\description{
  Computes the overlap between a pair of
  \code{GenomicRanges} objects
}
\details{
  Compute the overlap between a pair of
  \code{GenomicRanges} objects
}
\author{
  Giovanni d'Ario
}

