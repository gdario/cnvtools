\name{mapGeneToProbes}
\alias{mapGeneToProbes}
\title{mapGeneToProbes}
\usage{
  mapGeneToProbes(entrezid,
    markerRangesFile = paste("..", "compare_ge_and_cnv", "RData", "MarkerRanges.RData", sep = .Platform$file.sep))
}
\arguments{
  \item{entrezid}{character vector containing the ENTREZ
  IDs}
}
\value{
  a character vector containing the CNV probes.
}
\description{
  Map ENTREZ IDs to CNV probes
}
\details{
  Given a list of ENTREZ IDs, map them to the corresponding
  CNV probes.
}
\author{
  Giovanni d'Ario
}

