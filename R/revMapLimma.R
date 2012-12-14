##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title 
##' @param formula 
##' @param RevMappedCN 
##' @param clinicalInfo 
##' @param outputPrefix 
##' @param saveOutput 
##' @param match 
##' @param matchCol 
##' @return 
##' @author ubuntu
revMapLimma <- function(formula=NULL,
                        RevMappedCN=NULL,
                        clinicalInfo=commonAnnotationGEandCNV,
                        outputPrefix="out_",
                        saveOutput=TRUE,
                        match=TRUE,
                        matchCol="CNVsampleName2") {

    if(!require(limma))
	stop("Cannot find the 'limma' package")
    
    ## If necessary, match the colnames of the RevMappedCN matrix
    ## with the appropriate column of the clinical info
    if(match) {
        idx <- match(colnames(RevMappedCN), clinicalInfo[[matchCol]])
        if(any(is.na(idx)))
           stop(paste("ERROR: some patients do not match between",
                      "the column names of the RevMappedCN matrix",
                      "and the", matchCol, "column of the clinicalInfo"))
        clinicalInfo <- clinicalInfo[idx, , drop = TRUE]
        rownames(clinicalInfo) <- clinicalInfo[[matchCol]]
    }
        
    Data <- model.frame(formula = formula,
			data = clinicalInfo)
    design <- model.matrix(formula, data = Data)

    ## Match the copy number data and the design matrix rownames
    idx <- match(rownames(design), colnames(RevMappedCN))
    RevMappedCN <- RevMappedCN[, idx]

    ## By default consider all the chromosomes
    ## This might change in the feature
    
    revMapLimmaFit <- lmFit(RevMappedCN, design = design, )
    revMapLimmaFit <- eBayes(revMapLimmaFit)
    
    ## Make TopTable an object of class revMapLimma
    attr(revMapLimmaFit, "class") <- c("revMapLimmaFit", "list")
    attr(revMapLimmaFit, "formula") <- formula
    
    if(saveOutput) {
        fileName <- paste(outputPrefix, ".RData", sep = "")
        save(revMapLimmaFit, file = fileName)
    }
    
    return(invisible(revMapLimmaFit))
}


