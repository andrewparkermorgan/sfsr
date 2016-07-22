## io.R
## read and write SFS from disk

#' Read 'flattened' site frequency spectrum (SFS) from a file
#' 
#' @param ff filename
#' @param dims sample sizes (number of chromosomes) along each dimension of SFS
#' @param dtype what sort of frequences to expect; best let the procedure auto-detect
#' @param bootstraps logical; if \code{TRUE}, read bootstrap replicates (one per line) if they are present
#' @param repolarize logical; if \code{TRUE}, swap ancestral and derived states (see Details)
#' @param ... ignored
#' 
#' @return a site frequency spectrum (SFS): a k-dimensional array representing the joint frequencies of derived alleles in each of k populations
#' 
#' @details The SFS is expected to be provided in a text file as space-separated numbers (integers or floating-point), in row-major order.  
#' An \emph{unfolded} SFS (ie. polarized against the derived allele) is expected, although this may be relaxed in future
#'   
#' NB: Recent versions of \code{ANGSD} apparently get the ancestral and derived alleles backwards.  Use \code{repolarize = TRUE} to correct
#' this issue at runtime.
#' 
#' @export
read_sfs <- function(ff, dims, dtype = double, bootstraps = TRUE, repolarize = FALSE, ...) {
	
	#hack
	angsd <- FALSE
	
	## read comment lines with metadata
	message("Reading SFS metadata...")
	infile <- file(ff)
	open(infile)
	meta <- vector("list")
	mn <- character()
	j <- 1
	while (length(line <- readLines(infile, n = 1, warn = FALSE))) {
		
		if (!grepl("^#", line)) {
			break
		}
		
		vals <- stringr::str_match(line, "#+(\\w+)=(.+)")[1,2:3]
		vv <- unlist(stringr::str_split(vals[2], ","))
		mn <- c(mn, vals[1])
		vv.num <- suppressWarnings(as.numeric(vv))
		vv.bool <- suppressWarnings(as.logical(vv))
		if (!any(is.na(vv.num)))
			meta[[j]] <- vv.num
		else if (!any(is.na(vv.bool)))
			meta[[j]] <- vv.bool
		else
			meta[[j]] <- vv
		
		j <- j+1
		
	}
	close(infile)
	names(meta) <- mn
	
	if ("pops" %in% names(meta))
		names(dims) <- meta[["pops"]]
	
	## convenience wrapper for base::scan()
	read_one <- function(ff, ...) {
		line <- scan(ff, what = dtype(), comment.char = "#", quiet = TRUE,
					 nlines = 1, strip.white = TRUE, ...)
		return(line)
	}
	
	## now read spectrum itself
	message("Reading SFS: expecting ", prod(dims+1), " entries...")
	boots <- list()
	infile <- file(ff, "r")
	ii <- 1
	if (!angsd) {
		while(length(x <- read_one(infile)) || ii == 1) {
			if (length(x) == 0)
				next
			else if (ii > 1 && !bootstraps)
				break
			sfs <- aperm( array(x, dim = rev(dims+1)) )
			dimnames(sfs) <- lapply(dims, function(z) c(0,seq_len(z)))
			boots[[ii]] <- sfs
			ii <- ii+1
		}
	}
	else {
		while(length(x <- read_one(infile))) {
			if (ii > 1 && !bootstraps)
				break
			sfs <- array(x, dim = dims+1)
			dimnames(sfs) <- lapply(dims, function(z) c(0,seq_len(z)))
			boots[[ii]] <- sfs
			ii <- ii+1
		}
	}
	close(infile)
	
	## add metadata
	rez <- boots[[1]]
	message("Read ", length(rez), " items.")
	if (length(boots) > 1) {
		message("\t[ ", length(boots), " bootstrap replicates ]")
		attr(rez, "bootstraps") <- boots
	}
	if (length(meta)) {
		for (m in names(meta))
			attr(rez, m) <- meta[[m]]
	}
	
	if (repolarize)
		rez <- repolarize_sfs(rez)
	
	message("First SFS has ", sum(rez), " sites.")
	if (!inherits(rez, "sfs"))
		class(rez) <- c("sfs", class(rez))
	
	
	return(rez)
	
}

#' Convert SFS to vector representation for saving in a file
#' 
#' @param sfs an SFS
#' 
#' @return a vector containing the 'flattened' representation of the SFS, in row-major order
#' 
#' @export
flatten_sfs <- function(sfs) {
	
	## because numpy flattens row-wise and R column-wise
	as.vector(aperm(sfs))
	
}