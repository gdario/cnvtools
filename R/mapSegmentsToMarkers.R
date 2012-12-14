##' Finds the markers falling into a set of segments.
##'
##' Given a matrix of segments (the 'output' component of a DNAcopy
##' object) 
##' and a \code{GRanges} object containing the (unitary) ranges of the
##' markers on the MIP platform, the function returns the
##' \code{seg.means} (the average value of the segments) mapped to
##' the probes falling within the segment itself.
##' @title Finds the markers falling into a set of segments.
##' @param segs the 'output' component of a DNAcopy object containing
##' the smoothed and segmented values of either the copy number or the
##' allele ratio.
##' @param assayRanges a GRanges object containing the coordinates of
##' the probes. Such coordinates are represented as segments where the
##' start and the end are the same.
##' @param assayData a data frame containing the probe-level information 
##' @return a vector of segmentation mean values whose names are the
##' names of the probes falling within that segment.
##' @author Giovanni d'Ario
seg2assay <- function(segs, assayRanges, assayData) {
	
	## Sorting never hurts...
	segs <- segs[order(segs$loc.start, decreasing = FALSE), ]
	
	## Create a GRanges object containing the ranges of the segments
	SegRanges <- GRanges(seqnames = segs$chrom,
			ranges = IRanges(start = segs$loc.start,
					end = segs$loc.end))
	
	## Find the segments each marker falls in
	overlaps <- findOverlaps(assayRanges, SegRanges)
	
	## TODO: this error message is farly obscure
	if(queryLength(overlaps) != nrow(assayData))
		stop("Not all the assays have been mapped on the patient")
	
	overlaps <- as.matrix(overlaps)
	probes <- assayRanges@elementMetadata[[1]]
	out <- segs$seg.mean[overlaps[, 2]]
	names(out) <- probes[overlaps[, 1]]
	return(out)
}

##' Turns a mapping list into a matrix.
##'
##' The \code{mapSegmentsToMarkers} functions creates a list of segments
##' as long as are the patients in the dataset. The function then
##' operates on the list applying the \code{seg2assay} function and
##' returns another list. This function takes care of transforming that
##' list into a matrix.
##' @title Turns a mapping list into a matrix.
##' @param mapList a list object, output of the
##' \code{mapSegmentsToMarkers} function.
##' @return a matrix with as many columns as are thee patients and as
##' many rows as are the markers.
##' @author Giovanni d'Ario
revMapToMatrix <- function(mapList) {
	
	probes <- lapply(mapList, names)    # Extract probe names
	nProbes <- sapply(mapList, length)  # Extract number of probes
	
	## Which patient has the largest number of mapped probes?
	idxMax <- which.max(nProbes)
	maxProbes <- probes[[idxMax]]
	
	tmp <- lapply(mapList, function(x) {
				idx <- match(maxProbes, names(x))
				return(x[idx])
			})
	out <- do.call("cbind", tmp)
	
	## The 'out' matrix now contains some NAs, corresponding to the
        ### patients for which not all probes could be mapped.
        ## We remove the NAs
	idxNA <- apply(out, 1, function(x) any(is.na(x)))
	out <- out[!idxNA, ]
	return(out)
}

##' Reverse mapping of the segment values to the MIP markers.
##'
##' This function takes an object of class \code{CNVset} or
##' \code{DNAcopy} and a chromosome. For each patient, it assigns the
##' segmented copy number or allele ratio estimates to the probes that
##' fall inside the segment. This provides an "ultasmoothed" version
##' of the CN or AR estimates, with the advantage that now each sample
##' will have the same number of elements, thus allowing direct
##' comparisons by statistical tests and/or clustering.
##' @title Reverse mapping of the segment values to the MIP markers.
##' @param x an object of class \code{DNAcopy}.
##' @param chrom integer, the chromosome of interest (1-22)
##' @param cutoff integer representing the cutoff value for 'small'
##' copy number values. Default is 0.1. 
##' @return A matrix containing as many columns as are the samples
##' in the dataset and as many rows as are the markers mapping the
##' chromosome of interest.
##' @author Giovanni d'Ario
##' @export 
mapSegmentsToMarkers <- function(x, chrom=NULL) {
    
    if(inherits(x, "DNAcopy")) {
        ## Data <- x$data
        Segments <- x$output
    } else {
        stop("The x object must be of 'DNAcopy' class")
    }
    
    ## Subset the Data and the Segments to the chromosome of interest 
    ## if(!chrom %in% Data$chrom)
    ##    stop("Invalid chromosome name")
    
    ## Data <- Data[Data$chrom == chrom, ]
    Segments <- Segments[Segments$chrom == chrom, ]
        
    ## Subset the mipMarkers to the chromosome of interest
    mipMarkers$Chrom <- sub("chr0{0,1}", "", mipMarkers$Chrom)
    Markers <- mipMarkers[mipMarkers$Chrom == chrom, ]
    Markers <- Markers[order(Markers[["Chrom Pos (bp)"]],
                             decreasing = FALSE), ]
    
    ## Create a GRanges object containing the ranges of the Markers
    AssayRanges <- GRanges(seqnames = Markers$Chrom,
                           ranges = IRanges(start = Markers[["Chrom Pos (bp)"]],
                               end = Markers[["Chrom Pos (bp)"]]),
                           probeID = Markers[["Assay Name"]])
    
    ## Create a list of Segments, one component for each patient
    SegList <- split(Segments, Segments$ID)
    
    message(paste("I am reverse-mapping the segments for chromosome",
                  chrom))
    
    out <- lapply(SegList,
                  seg2assay,
                  assayRanges = AssayRanges,
                  assayData = Markers)
    
    ## The revMapToMatrix function takes care of selecting the common
    ## markers and put everything in a matrix format
    out <- revMapToMatrix(out)
    
    return(out)
}

##' Check that the probes in a matrix of reverse-mapped segments are
##' in the same order as in the genome.
##'
##' In principle all the reverse-mapped segments should give rise to
##' matrices
##' that are ordered in increasing order of the genomic location, but
##' some manipulations could change this order. \code{testOrder} checks
##' that the order is preserved. IMPORTANT: I noticed that some probes 
##' are missing in the rev-maps. This is strange and should be further
##' investigated.
##' @title testOrder
##' @param x a list of matrices of rev-mapped segments.
##' @param chrom numeric or character, indicating the chromosome of
##' interest.
##' X and Y chromosomes are excluded from the analysis
##' If NULL the default file (orderedMipAssays.RData) is loaded
##' @return a matrix containing rev-mapped segments ordered by genomic
##' location
##' @author Giovanni d'Ario
testOrder <- function(x=NULL, chrom=NULL) {
	
	## is x an object of class 'list'?
	if(!inherits(x, "list"))
		stop("Error: x must be a list (one component per chromosome)")
	
	## Chrom can be either a character or a number. If it is a
	## character it must be convertible into a number simply
	## by applying 'as.numeric'
	if(is.character(chrom))
		chrom <- as.numeric(chrom)
	if(is.na(chrom))
		stop("Error: invalid chromosome entered")
	
	## Extract the rev-mapped segments
	RevMap <- x[[chrom]]
		
	## remove che 'chr0' from mipMarkers$Chromosome
	mipMarkers$Chromosome <- sub("^chr0{0,1}", "",
			mipMarkers$Chromosome)
	mipMarkers$Chromosome <- as.numeric(mipMarkers$Chromosome)
	## Remove the NAs due to the X and Y chromosomes
	mipMarkers <- mipMarkers[!is.na(mipMarkers$Chromosome), ]
	
	## Subset the mipMarkers to the chrom of interes
	markersInChrom <- subset(mipMarkers, Chromosome == chrom)
	
	## FIXME: some markers are missing in the rev-mapped matrix
	## WHY?
	idx <- markersInChrom$Assay_Name %in% rownames(RevMap)
	markersInChrom <- markersInChrom[idx, ]
	
	## Make sure they are ordered by position
	if(!all.equal(markersInChrom$Chrom_Position,
			sort(markersInChrom$Chrom_Position)))
		markersInChrom <- markersInChrom[order(markersInChrom$Chrom_Position),]
	
	## Check that the rev-mapped segments appear in the same order as
	## the sorted markers
	if(!all.equal(rownames(RevMap), markersInChrom$Assay_Name)) {
		message("The probes in 'x' were not ordered. Fixing it")
		RevMap <- RevMap[match(markersInChrom$Assay_Name,
						rownames(RevMap)), ]
	}
	return(RevMap)
}

##' Create a list of allele ratios matrices with 22 components (the
##' number of autosomal chromosomes), ordered as the corresponding
##' copy number reverse-mapped matrices.
##'
##' Copy number estimates have been segmented and then the segments have
##' been reverse-mapped to matricexs in order to make the use of heatmaps
##' and ordinary statistical tools available. I have however noticed that
##' segmentation has a possibly detrimental effect when applied to the
##' allele ratio estimates. This function takes the CNA object
##' @title makeMatchedListOfARs
##' @param ARobject a CNA object (tipically the 'DNAcopyAlleleRatios'
##' object)
##' containing the smoothed but *not* the segmented allele ratios for
##' each sample and each probe
##' @param CNobject a list of matrices of reverse-mapped CN segments.
##' @return a list of matrices containing the allele ratios in the same
##' row and column order as the corresponding CN rev-mapped list.
##' @author Giovanni d'Ario
makeMatchedListOfARs <- function(ARobject=NULL, CNobject=NULL) {
	
	Lst <- list()
	tmp <- ARobject[, -c(1:3)]
	rownames(tmp) <- ARobject$assay_name
	
	for(i in seq_along(CNobject)) {
		## Make sure the patients appear in the same order in the
		## CN and the AR data
		CN <- CNobject[[i]]
		
		tmp <- tmp[, match(colnames(CN), colnames(tmp))]
		idx <- match(rownames(CN), rownames(tmp))
		AR <- tmp[idx,]
		Lst[[i]] <- as.matrix(AR)
	}
	
	return(Lst)
}

