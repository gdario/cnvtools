##################################################################
getChromAndCoords <- function (probes) {

    if(!exists("MipMarkers"))
        stop("I cannot find the MIP probe-level info (MipMarkers)")
    if (!inherits(probes, "character")) 
        stop("Error: 'probes' must be a character vector")

    idx <- match(probes, MipMarkers$probe_id)
    coords <- MipMarkers$chrom_pos[idx]
    chrom <- MipMarkers$chrom[idx]

    M <- cbind(chrom, coords)
    colnames(M) <- c("chrom", "coords")
    rownames(M) <- probes
    return(M)
}

##################################################################
getCoefsAndPvals <- function(revMapLimmaFit=NULL,
                             p.adj=NULL) {
    
    f <- attr(revMapLimmaFit, "formula") 
    hasIntercept <- !as.logical(length(grep("- *1", f)))
    
    if(hasIntercept) {
        idx <- match("(Intercept)",
                     colnames(revMapLimmaFit$coefficients))
        Coef <- revMapLimmaFit$coefficients[,-idx, drop = FALSE]
        Pv <- revMapLimmaFit$p.value[, -idx, drop = FALSE]
    } else {
        Coef <- revMapLimmaFit$coefficients
        Pv <- revMapLimmaFit$p.value
    }
    
    ## Make sure that Coef and Pv have the probes in the
    ## same order (it should always be the case).
    if(!all.equal(rownames(Coef), rownames(Pv))) {
        idx <- match(rownames(Coef), rownames(Pv))
        Pv <- Pv[idx, ]
    }

    ## Get the chromosome and the genomic coordinates
    rn <- rownames(Coef)
    chromAndCoords <- getChromAndCoords(rn)

    ## Make sure the the Coef, Pv and chromAndCoords rows are
    ## in the same order.
    if(!all.equal(rownames(Coef), rownames(chromAndCoords))) {
        idx <- match(rownames(Coef), rownames(chromAndCoords))
        chromAndCoords <- chromAndCoords[idx, ]
    }

    ## Put the chromosome and the genomic coordinates as the
    ## first two columns, for easier elimination in the plots
    Coef <- cbind(chromAndCoords, Coef)
    Pv <- cbind(chromAndCoords, Pv)

    ## Make sure that the chromosomes and the genomic coordinates
    ## are in increasing order
    idx <- order(chromAndCoords[,1], chromAndCoords[,2])
    if(!all.equal(idx, 1:length(idx))) {
        Coef <- Coef[idx, ]
        Pv <- Pv[idx, ]
    }

    ## Turn the matrices into data frames. This is
    ## particularly important for adjustPv, which is
    ## supposed to work on data.frames, and will return
    ## the wrong result (after a long time), if 
    ## provided with a matrix.
    Pv <- as.data.frame(Pv)
    Coef <- as.data.frame(Coef)
    
    ## Adjust the p-values
    AdjPv <- adjustPvals(Pv, p.adj=p.adj)
    
    ## Create the output list and return it
    coefsAndPvals <- list(coefficients = Coef,
                pvalues = AdjPv)
    return(coefsAndPvals)
}

##################################################################
adjustPvals <- function(Pv=NULL, p.adj=NULL) {
    idx <- match(c("chrom", "coords"), colnames(Pv))
    if(p.adj != "none") {
        message(paste("Adjusting the p-values using the",
                      p.adj, "method. Please wait."))
        PvMat <- Pv[-idx]
        tmp <- sapply(PvMat, p.adjust, method = p.adj)
        Pv[-idx] <- tmp
    }
    return(Pv)
}

##################################################################
plotPvalues <- function(CoefAndPvals=NULL,
                        zoom=NULL,
                        pv.lim=NULL,
                        legend=NULL) {
    
    if(!is.null(zoom))
        CoefAndPvals <- zoomIn(CoefAndPvals)

    Pvals <- CoefAndPvals$pvalues
    
    ## The number of rows of the TopTable object is
    ## used instead of the genomic coordinates.
    nr <- nrow(Pvals)
    idx <- match(c("chrom", "coords"), colnames(Pvals))
    
    ## Negative log10 of the p-values
    Pvals[-idx] <- -log10(Pvals[-idx])
    rangeLpv <- sapply(Pvals[-idx], range)

    ## RangeLpv is a matrix, we take the range of the ranges
    ## to obtain a vector
    rangeLpv <- range(rangeLpv)
    
    ## Plot the empty plot
    plot(1:nr, ylim = rangeLpv, type = "n",
         axes = FALSE, xlab = NA, ylab = NA,
         main = "-log10(P-value)")

    ## Add the chromosome rectangles and the axes
    plotRectanglesAndAxes(Pvals, nr,
                          xtext = "Chromosome",
                          ytext = "-log10(p-value)",
                          ybottom = rangeLpv[1],
                          ytop = rangeLpv[2])
    
    ## Plot the -log10(p-values)
    PvalsMat <- Pvals[, -idx, drop=FALSE]
    for(i in 1:ncol(PvalsMat)) {
        lines(1:nr, PvalsMat[,i],
                    type = "s",
                    col = i)
    }
    
    ## Add a line showing the limit p-value
    segments(x0 = 1, x1 = nr,
             y0 = -log10(pv.lim),
             y1 = -log10(pv.lim),
             lty = 2, col = "blue")

    addLegend(Pvals, legend)
}

##################################################################
plotCoefs <- function(CoefAndPvals=NULL,
                      zoom=NULL,
                      legend=NULL) {

    if(!is.null(zoom))
        CoefAndPvals <- zoomIn(CoefAndPvals)

    Coefs <- CoefAndPvals$coefficients
    
    ## The number of rows of the TopTable object is
    ## used instead of the genomic coordinates.
    nr <- nrow(Coefs)
    idx <- match(c("chrom", "coords"), colnames(Coefs))
    
    ## Get the range of the coefficients
    rangeCoef <- sapply(Coefs[-idx], range)

    ## RangeLpv is a matrix, we take the range of the ranges
    ## to obtain a vector
    rangeCoef <- range(rangeCoef)
   
    ## Empty plot
    plot(1:nr, ylim = rangeCoef, type = "n",
         axes = FALSE, xlab = NA, ylab = NA,
         main = "Estimated coefficients")

    ## Plot the rectangles and the axes
    plotRectanglesAndAxes(Coefs, nr,
                          xtext = "Chromosome",
                          ytext = "Estimated coefficients",
                          ybottom = rangeCoef[1],
                          ytop = rangeCoef[2])

    ## Plot the estimated coefficients
    CoefMat <- Coefs[, -idx, drop=FALSE]
    for(i in 1:ncol(CoefMat)) {
        lines(1:nr, CoefMat[,i],
              type = "s",
              col = i)
    }

    ## Add the 0 copy number line
    segments(x0 = 1, x1 = nr,
             y0 = 0,
             y1 = 0,
             lty = 2, col = "blue")    

    ## Add the legends
    addLegend(Coefs, legend)
}

##################################################################
getRectangleCoords <- function(obj, nr) {

    ## Create the chromosome rectangles
    dd <- c(0, diff(obj$chrom))
    xbk <- which(dd == 1)
    xleft <- c(0, xbk)
    xright <- c(xbk, nr)
    
    ## Position of the chromosome indicators
    labPos <- .5 * (xleft + xright)

    ## Output list containing the coordinates
    rectCoords <- list(dd = dd,
                       xbk = xbk,
                       xleft = xleft,
                       xright = xright,
                       labPos = labPos)
    return(rectCoords)
}

##################################################################
plotRectanglesAndAxes <- function(obj,
                                  nr,
                                  xtext,
                                  ytext,                                  
                                  ybottom,
                                  ytop) {
    
    ## Compute the rectangle coordinates
    rectCoords <- getRectangleCoords(obj, nr=nr)

    ## Add the chromosome rectangles
    rect(xleft = rectCoords$xleft,
         xright = rectCoords$xright,
         ybottom = ybottom,
         ytop = ytop,
         col = c("white", "lightGrey"))

    ## Set the axes
    setAxes(rectCoords, xtext, ytext)
}

##################################################################
addLegend <- function(obj, legend) {
    idx <- match(c("chrom", "coords"), colnames(obj))
    objMat <- obj[, -idx, drop = FALSE]
    if(is.null(legend))
        legend <- colnames(objMat)
    legend(x = "topleft",
           bty = "n",
           legend = legend,
           col = 1:ncol(objMat),
           text.col=1:ncol(objMat))
}

##################################################################
setAxes <- function(rectCoords=NULL,
                    xtext=NULL,
                    ytext=NULL) {
    ## Set the axes
    axis(1, at = rectCoords$labPos, labels = 1:22)
    mtext(xtext, side = 1, line = 3)
    axis(2)
    mtext(ytext, side = 2, line = 3)
}

##################################################################
zoomIn <- function(CoefAndPvals, zoom) {
    idx <- with(CoefAndPvals$coefficients,
                chrom == zoom[1] &
                coords >= zoom[2] &
                coords <= zoom[3])
    CoefAndPvals <- lapply(CoefAndPvals, function(x)
                           return(x[idx, ]))
    return(CoefAndPvals)
}

##################################################################
plot.revMapLimmaFit <- function(revMapLimmaFit=NULL,
                                zoom=NULL,
                                pv.lim=0.05,
                                legend=NULL,
                                p.adj = c("none",
                                    "holm",
                                    "hochberg",
                                    "hommel",
                                    "bonferroni",
                                    "BH",
                                    "BY",
                                    "fdr")) {
    
    p.adj <- match.arg(p.adj)

    ## Extract the coefficients and the (possibly) adjusted
    ## p-values
    CoefAndPvals <- getCoefsAndPvals(revMapLimmaFit,
                                     p.adj = p.adj)
    op <- par(mfrow = c(2,1))

    ## Plot the p-values
    plotPvalues(CoefAndPvals = CoefAndPvals,
                zoom = zoom,
                legend=legend,
                pv.lim = pv.lim)

    ## Plot the coefficient estimates
    plotCoefs(CoefAndPvals = CoefAndPvals,
              legend = legend,
              zoom = zoom)

    on.exit(par(op))
}

##################################################################
## plot.revMapLimma <- function(revMapLimmaFit,
##                              zoom=NULL,
##                              p.adj = c("none",
##                                  "holm",
##                                  "hochberg",
##                                  "hommel",
##                                  "bonferroni",
##                                  "BH",
##                                  "BY",
##                                  "fdr"),
##                              pv.lim=0.05) {

##     ## Select the type of correction
##     p.adj <- match.arg(p.adj)

##     ## Get the coefficients and the p-values
##     CoefAndPvals <- getCoefsAndPvals(revMapLimmaFit = revMapLimmaFit,
##                                      p.adj = p.adj)

##     plotPvAndCoef(CoefAndPvals,
##                   zoom = zoom,
##                   which.pv = pv,
##                   pv.lim = pv.lim)
## }
