\name{mapEntrezIdToSegments}
\alias{mapEntrezIdToSegments}
\title{mapEntrezIdToSegments}
\usage{
  mapEntrezIdToSegments(entrezid,
    segmentFile = "../common_cnv_and_ge_annotation/RData/cnv_segments_304_patients.RData")
}
\arguments{
  \item{entrezid}{character vector containing the ENTREZ
  IDs}

  \item{segmentFile}{string containing the full path to a
  file containing the output of a segmentation algorithm,
  like DNAcopy (the output slot).}
}
\value{
  a data frame containing the ENTREZ IDs and the segments
  that overlap with the genomic coordinates of the genes.
}
\description{
  Map a set of ENTREZ IDs to the segments contained in a
  'segment file'.
}
\details{
  Given a set of ENTREZ IDs, find the segments in a
  segmented data set overlapping with the corresponding
  genes.
}
\author{
  Giovanni d'Ario
}

