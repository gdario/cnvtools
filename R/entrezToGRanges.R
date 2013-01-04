##' Converts a set of ENTREZ IDs into genomic ranges (HG18)
##'
##' Converts a set of ENTREZ IDs into genomic ranges (HG18)
##' @title entrezToGRanges
##' @param entrezid character vector containing the ENTREZ IDs of the
##' genes of interest.
##' 
##' @return and object of class \code{GenomicRanges}.
##' @author Giovanni d'Ario
##' @export
entrezToGRanges <- function(entrezid) {
    require(biomaRt)
    require(GenomicRanges)
    
    ## Get the ensembl corresponding to the hg18
    ensembl <- useMart("ENSEMBL_MART_ENSEMBL",
                       dataset="hsapiens_gene_ensembl",
                       host="may2009.archive.ensembl.org",
                       path="/biomart/martservice",

                       archive=FALSE)
    
    eid <- unique(entrezid)
    GeneCoords <- getBM(filters = "entrezgene",
                        attributes = c("chromosome_name",
                            "start_position",
                            "end_position",
                            "entrezgene"),
                        values = eid,
                        mart = ensembl)
    
    GeneRanges <- GRanges(seqname = GeneCoords$chromosome_name,
                      ranges = IRanges(
                          GeneCoords$start_position,
                          end = GeneCoords$end_position),
                      ## strand = ifelse(GeneCoords$strand == 1,
                      ##     "+", "-"),
                      entrezgene = as.character(
                          GeneCoords$entrezgene))
    return(GeneRanges)
}
