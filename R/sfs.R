## sfs.R
## methods for reading, manipulating and testing site frequency spectra (SFS)

## check if an SFS has bootstrap replicates attached
.hasboots <- function(sfs) !is.null(attr(sfs, "bootstraps"))

## apply function over the bootstrap replicates attached to an SFS
.sapply.boot <- function(sfs, f, ...) sapply(attr(sfs, "bootstraps"), f, ...) 
.lapply.boot <- function(sfs, f, ...) lapply(attr(sfs, "bootstraps"), f, ...)

## check that all SFS have compatible sets of bootstraps
.compatible.boots <- function(sfs) {
	if (all(sapply(sfs, .hasboots))) {
		lens <- sapply(sfs, function(f) length(attr(f, "bootstraps")))
		return(length(unique(lens)) == 1)
	}
	else {
		return(FALSE)
	}
}

## get a specific bootstrap replicate
.get.boot <- function(sfs, ii) {
	attr(sfs, "bootstraps")[[ii]]
}

#' Discard bootstrap replicates from an SFS
#' 
#' @param sfs an SFS
#' @return an SFS of same dimensions with any bootstrap replicates discarded
#' 
#' @export
noboots <- function(sfs) { 
	attr(sfs, "bootstraps") <- NULL
	return(sfs)
}

#' Get an SFS with reduced dimensions
#' 
#' @param sfs an SFS
#' @param dims indices over which to sum up allele frequencies
#' @param persite logical; if \code{TRUE}, divide by the total number of sites
#' @param ... ignored
#' 
#' @return an SFS of reduced dimension, produced by "marginalizing out" (ie. just summing over) one or more populations
#' 
#' @export
marginalize_sfs <- function(sfs, dims = 1, persite = FALSE, ...) {
	
	.marginalize <- function(x) {
		xx <- apply(x, dims, sum)
		if (persite)
			xx <- xx/sum(xx)
		if (length(dim(xx)))
			dimnames(xx) <- dimnames(sfs)[dims]
		else
			names(xx) <- dimnames(sfs)[[dims]]
		return(xx)
	}
	
	rez <- .marginalize(sfs)
	if (.hasboots(sfs)) {
		attr(rez, "bootstraps") <- .lapply.boot(sfs, .marginalize)
	}
	class(rez) <- c("sfs", class(rez))
	return(rez)
}

#' Aggreage an n-dimensional SFS into a single dimension
#' 
#' @param sfs an SFS
#' @param ... ignored
#' 
#' @return an one-dimensional SFS
#' 
#' @export
aggregate_sfs <- function(sfs, ...) {
	
	.aggregate <- function(x) {
		if (ndim_sfs(x) == 1) {
			return(x)
		}
		else {
			ii <- which(!is.na(x), arr.ind = TRUE)
			jj <- rowSums(ii-1)
			new.sfs <- rep(0, sum(dim_sfs(x))+1)
			for (k in seq_along(jj)) {
				new.sfs[ jj[k]+1 ] <- new.sfs[ jj[k]+1 ] + x[ ii[k,, drop = FALSE] ]
			}
			return(new.sfs)
		}
	}
	
	rez <- .aggregate(sfs)
	names(rez) <- seq(0, length(rez)-1)
	if (.hasboots(sfs))
		.lapply.boot(sfs, .aggregate)
	
	if (!inherits(rez, "sfs"))
		class(rez) <- c("sfs", class(rez))
	
	return(rez)
	
}

print.sfs <- function(sfs, ...) {
	cat("Site frequency spectrum with dimensions [", paste(dim_sfs(sfs), collapse = " "), "]\n")
	cat("\t", round(sum(sfs)), "sites\n")
	if (.hasboots(sfs))
		cat("\t", length(attr(sfs, "bootstraps")), "bootstrap replicates\n")
	cat("\n")
	attr(sfs, "bootstraps") <- NULL
	sfs <- round(sfs)
	NextMethod("print")
}

## project SFS from one sample size to another
## logic borrowed from Dadi: <https://github.com/paulirish/dadi/blob/master/dadi/Spectrum_mod.py>
project_sfs <- function(sfs, newdims, ...) {
	
	.proj.coef <- function(nto, nfrom, nhit) {
		phit <- seq(0, nto)
		exp( lchoose(nto, nhit) + lchoose(nfrom-nto, nhit-phit) - lchoose(nfrom, nhit) )
	}
	
	## the projected SFS, must have smaller size than parent
	rez <- array(0, dim = newdims)
	stopifnot(all(newdims <= dim(sfs)))
	
	## project along one axis
	.project.one <- function(n, axis) {
		
		from <- dim(sfs)[axis]
		for (hits in seq(1, from)) {
			least <- max(n-(from-hits), 0)
			most <- min(hits, n)
			print(c(least,most))
			for (j in seq(least, most)) {
				idx <- c(1,1,1)
				idx[axis] <- hits
				## UGGH
			}
		}
	}
	
	for (a in seq_along(newdims)) {
		.project.one(newdims[a], a)
	}
	
	return(rez)
	
}

## Force an SFS to be a matrix
matrify <- function(x, ...) {
	ndim <- length(dim(x))
	if (!(ndim > 1)) {
		x <- matrix(x, ncol = length(x))
	}
	if (!inherits(x, "sfs"))
		class(x) <- c("sfs", class(x))
	return(x)
}

ndim_sfs <- function(x, ...) {
	
	if (!length(dim(x)))
		return(1)
	else if (min(dim(x)) == 1)
		return(1)
	else
		return(length(dim(x)))
	
}

dim_sfs <- function(x, ...) {
	dims <- dim(x)
	if (!length(dims))
		return(length(x)-1)
	else if (min(dims) == 1)
		return(dims[ dims > 1 ]-1)
	else
		return(dims-1)
}

## check whether a set of SFS have all same dimensions
same_size <- function(...) {
	x <- list(...)
	if (is.list(x[[1]]))
		x <- x[[1]]
	dims <- lapply(x, dim_sfs)
	return( length(unique(dims)) == 1 )
}

#' Zero out frequency classes of fixed for ancestral or derived allele, leaving only polymorphic classes
#' 
#' @param sfs an SFS
#' @param ... ignored
#' 
#' @return the original SFS with fixed bins (frequency zero or N) set to zero
#' 
#' @export
mask_corners <- function(sfs, ...) {
	
	sfs <- matrify(sfs)
	
	ndims <- length(dim(sfs))
	sfs[ t(rep(1, ndims)) ] <- 0
	sfs[ t(dim(sfs)) ] <- 0
	attr(sfs, "masked") <- TRUE
	if (!inherits(sfs, "sfs"))
		class(sfs) <- c("sfs", class(sfs))
	return(sfs)
	
}

mask_singletons <- function(sfs, ...) {
	
	sfs[2] <- 0
	return(sfs)
	
}

#' Reverse the polarity (ancestral vs derived) of an SFS.
#' 
#' @param sfs an SFS
#' @param ... ignored
#' 
#' @details It appears that recent versions of the \code{ANGSD} software get the polarity of the SFS backwards,
#' 	at least in their output.  The distinction is not important for simple measures of diversity, but for tests
#' 	based on polymorphism vs divergence (McDonald-Kreitman, Hudson-Kreitman-Aguade) it matters.
#' 
#' @export
repolarize_sfs <- function(sfs, ...) {

	.repolarize <- function(x) {
		x <- matrify(x)
		ii <- which(!is.na(x), arr.ind = TRUE)
		xx <- array(x[ apply(ii, 2, rev) ], dim = dim(x))
		return(xx)
	}
	
	rez <- .repolarize(sfs)
	dimnames(rez) <- dimnames(sfs)
	if (.hasboots(sfs)) {
		attr(rez, "bootstraps") <- .lapply.boot(sfs, .repolarize)
	}
	
	class(rez) <- c("sfs", class(rez))
	return(rez)
	
}

## internal workhorse for calculating Watterson estimator
.theta.w <- function(sfs, persite = FALSE, ...) {
	
	n <- length(sfs)-1
	x <- mask_corners(t(as.matrix(sfs)))
	S <- sum(seq(1, n-1)^-1)
	theta <- sum(x)/S
	if (persite)
		return(theta/sum(sfs))
	else
		return(theta)
	
}

#' Calculate Watterson's estimator of theta
#' 
#' @return Watterson's \eqn{\theta}, the population-scaled mutation rate
#' 
#' @export
#' @rdname diversity_stats
theta_w <- function(sfs, persite = FALSE, ...) {
	
	if (length(dim(sfs)) > 1)
		stop("Watterson's estimator is only defined for a 1d SFS -- try marginalizing first.")
	
	if (.hasboots(sfs)) {
		.sapply.boot(sfs, .theta.w, persite = persite, ...)
	}
	else {
		.theta.w(sfs, persite, ...)
	}
}

## internal workhorse for computing avg pairwise differences
.pairwise.div <- function(sfs, persite = FALSE, ...) {
	denom <- sum(sfs)
	ndim <- ndim_sfs(sfs)
	n <- dim_sfs(sfs)
	if (ndim < 2) {
		# 1d sfs; the easy case
		S <- (0:n)*(n:0)
		pihat <- sum(sfs*S/choose(n,2))
	}
	else if (ndim == 2) {
		# 2d sfs; the harder case
		n <- dim_sfs(sfs)
		sfs <- mask_corners(sfs)
		S <- which(!is.na(sfs), arr.ind = TRUE)-1
		
		# d = a2*(n1-a1) + a1*(n2-a2)
		a1 <- outer(0:n[1], rep(1, n[2]+1))
		a2 <- outer(rep(1, n[1]+1), 0:n[2])
		ones <- outer(rep(1, n[1]+1), rep(1,n[2]+1))
		n1 <- n[1]*ones
		n2 <- n[2]*ones
		d <- sum(sfs*(a2*(n1-a1)+a1*(n2-a2)))
		pihat <- d/prod(n)
	}
	else {
		stop("No support for more than 2d SFS yet.")
	}
	if (persite)
		return(pihat/denom)
	else
		return(pihat)
}

#' Calculate the Tajima's pi-hat estimator of theta (avg pairwise differences)
#' 
#' @return value of \eqn{\hat{\pi} = \sum_{i>j}^N{d_{ij}}}, the average number of pairwise differences (aka expected heterozygosity)
#' 
#' @references
#' Tajima F (1983) Evolutionary relationship of DNA sequences in finite populations. Genetics 105: 437-360.
#' 
#' @export
#' @rdname diversity_stats
theta_pi <- function(sfs, persite = FALSE, ...) {
	
	if (ndim_sfs(sfs) > 1)
		stop("The theta_pi statistic is only defined for a 1d SFS -- try marginalizing first.")
	
	if (.hasboots(sfs)) {
		.sapply.boot(sfs, .pairwise.div, persite = persite, ...)
	}
	else {
		.pairwise.div(sfs, persite, ...)
	}
}

#' Calculate the Fu and Li's estimator of theta
#' 
#' @return value of \eqn{\hat{\theta_{\zeta}}}, the number of derived singletons in a sample
#' 
#' @references
#' Fu YX and Li WH (1993) Statistical tests of neutrality of mutations. Genetics 133: 693-709.
#' 
#' @export
#' @rdname diversity_stats
theta_zeta <- function(sfs, persite = FALSE, ...) {
	
	if (ndim_sfs(sfs) > 1)
		stop("The theta_zeta statistic is only defined for a 1d SFS -- try marginalizing first.")
	
	.tzeta <- function(x) unname(x[2])/ifelse(persite, sum(x), 1)
	if (.hasboots(sfs))
		.sapply.boot(sfs, .tzeta)
	else
		.tzeta(sfs)
	
}

#' Calculate the d_xy, the average pairwise divergence
#' 
#' @param sfs an SFS
#' @param persite logical; if \code{TRUE}, divide by number of sites
#' @param ... ignored
#' 
#' @return average pairwise divergence: if one population, to the outgroup (assuming unfolded, polarized SFS);
#' if two populations, between them
#' 
#' @export
#' @rdname diversity_stats
d_xy <- function(sfs, persite = FALSE, ...) {
	
	.dxy.1d <- function(x) sum((seq_along(x)-1)*x)/(length(x)-1)
	.dxy <- function(x) {
		denom <- if (!persite) 1 else sum(x)
		if (ndim_sfs(x) > 1)
			.pairwise.div(x, persite = persite)
		else
			.dxy.1d(x)/denom
	}
	
	sfs <- matrify(sfs)
	if (.hasboots(sfs))
		rez <- .sapply.boot(sfs, .dxy)
	else
		rez <- .dxy(sfs)
	
	return(rez)
	
}

#' Summarize SFS into 5 bins on each axis
#' 
#' @param sfs an SFS
#' @param ... ignored
#' 
#' @return an SFS of reduced dimension, which summarizes the full SFS into bins corresponding to singletons, doubletons, all other polymorhisms
#' 	(of higher sample frequency) and fixed differences.
#' 	
#' @details The "coarse-grained" SFS produced by this function is useful as a lower-dimensinoal summary statistic in ABC procedures.
#' 
#' @references
#' 	Tellier A, et al. (2011) Estimating parameters of speciation models based on refined summaries of the joint site-frequency spectrum.
#' 		PLoS One 6: e18155.
#' 
#' @export
coarsen_sfs <- function(sfs, ...) {
	
	idx.from <- as.matrix(expand.grid(lapply(dim(sfs), seq_len)))
	idx.to <- apply(idx.from, 2, function(f) { 
		f1 <- f
		maxsz1 <- pmin(4, max(f))
		maxsz2 <- pmin(5, max(f))
		interior <- f > 3 & f < max(f)
		margin <- f == max(f)
		f1[interior] <- maxsz1
		f1[margin] <- maxsz2
		return(f1)
	})
	
	#print(sum(sfs))
	#print(sum(sfs[ idx.from ]))
	#print(apply(idx.to, 2, max))
	
	rez <- array(0, pmin(5, dim(sfs)))
	for (row in seq_len(nrow(idx.to))) {
		rez[ idx.to[row,,drop = FALSE] ] <- rez[ idx.to[row,,drop = FALSE] ] + sfs[ idx.from[row,,drop = FALSE] ]
	}
	return(rez)

}

#' Count fixed versus polymorphic sites in each of K spectra
#' 
#' @param x an SFS, or list thereof
#' @param y another SFS, if \code{x} is not a list
#' @param bootstrap logical; if \code{TRUE}, make the contingency table for any bootstrap replicates which are present
#' @param persite logical; if \code{TRUE}, divide by total number of sites 
#' @param ... ignored
#' 
#' @return a 2 x K contingency table with count of fixed vs polymorphic sites in K spectra
#' 
#' @export
fixpoly <- function(x, y = NULL, persite = FALSE, bootstrap = FALSE, ...) {
	
	if (is.list(x)) {
		sfs <- x
		if (!is.null(y))
			sfs[[ length(sfs)+1 ]] <- y
	}
	else if (!is.null(y)) {
		sfs <- list(sites1 = x, sites2 = y)
	}
	else {
		sfs <- list(sites1 = x)
	}
	if (!same_size(sfs))
		stop("All SFS to be compared must have same dimensions.")
	
	rez <- sapply(sfs, .fixpoly, persite = persite)
	if (.compatible.boots(sfs) && bootstrap) {
		attr(rez, "bootstraps") <- lapply(seq_along(attr(sfs[[1]], "bootstraps")),
										  function(f) sapply(sfs, function(s) .fixpoly(.get.boot(s, ii = f), persite = persite)))
	}
	
	return(rez)
	
}

# .fixpoly <- function(x, ...) {
# 	
# 	x <- matrify(x)
# 
# 	idx <- which(!is.na(x), arr.ind = TRUE)
# 	maxd <- dim(x)
# 	poly <- apply(idx, 1, function(row) any(row < maxd) & any(row > 1))
# 	fixed <- apply(idx, 1, function(row) any(row == maxd & row > 1)) & !poly
# 	nfix <- sum(x[ idx[ fixed,, drop = FALSE] ])
# 	npoly <- sum(x[ idx[ poly,, drop = FALSE] ])
# 	
# 	return(c(fixed = nfix, polymorphic = npoly))
# 	
# }

.fixpoly <- function(x, persite = FALSE, ...) {
	
	nsites <- sum(x)
	if (ndim_sfs(x) == 1) {
		div <- unname(x[ length(x) ])
		poly <- unname(sum(x)-x[1]-div)
	}
	else {
		x <- mask_corners(matrify(x))
		div <- sum(x[,1]) + sum(x[1,])
		poly <- sum(x)-div
	}
	
	rez <- c(div = div, poly = poly)
	if (persite)
		rez <- rez/nsites
	return(rez)
	
}

## innards of the HKA test
.hka_test <- function(x, scale = 1.0, ...) {
	
	.Cn <- function(n) sum(1/seq_len(n-1))
	.Cnn <- function(n) sum(1/(seq_len(n-1)^2))
	
	## constants determined by sample size only
	ndim <- ndim_sfs(x[[1]])
	
	## case 1: 2 populations
	if (ndim > 1) {
		
		ns <- dim_sfs(x[[1]])
		na <- .Cn(ns[1])
		nb <- .Cn(ns[2])
		
		## get some preliminary quantities
		div <- sapply(x, d_xy)
		sa <- (sapply(x, function(f) .segsites(marginalize_sfs(f, 1))))
		Sa <- sum(sa)
		sb <- (sapply(x, function(f) .segsites(marginalize_sfs(f, 2))))
		Sb <- sum(sb)
		
		## solve least-squares eqns for coalescent pararameters
		fhat <- (Sb*na)/(Sa*nb)
		omega <- Sa/na
		Tbig <- sum(div)/omega - (1+fhat)/2
		theta <- (sa+sb+div)/(Tbig + (1+fhat)/2 + na + fhat*nb)
		
		## compute expectations
		esa <- theta*na
		esb <- theta*fhat*nb
		ed <- theta*(Tbig + (1+fhat)/2)
		vsa <- esa + (theta^2)*.Cnn(ns[1])
		vsb <- esb + ((fhat*theta)^2)*.Cnn(ns[2])
		vd <- ed + (theta*(1+fhat)/2)^2
		
		## test statistic, degrees of freedom
		Xsq <- sum(((sa-esa)^2)/vsa) + sum(((sb-esb)^2)/vsb) + sum(((div-ed)^2)/vd)
		df <- pmax(2*length(sa)-2, 1)
		
		## compute residuals
		obs <- cbind(sa, sb)
		expect <- cbind(esa, esb)
		resids <- (obs-expect)/sqrt(expect)
		
	}
	## case 2: 1 population + outgroup
	else if (ndim == 1) {
		
		ns <- dim_sfs(x[[1]])
		na <- .Cn(ns[1])
		
		## get some preliminary quantities
		div <- sapply(x, d_xy)
		sa <- sapply(x, .segsites)
		Sa <- sum(sa)
		
		## solve least-squares eqns for coalescent pararameters
		fhat <- NA
		omega <- Sa/na
		Tbig <- sum(div)/omega - 1
		theta <- (sa+div)/(Tbig + 1 + na)
		
		## compute expectations
		esa <- theta*na
		ed <- theta*(Tbig + 1)
		vsa <- esa + (theta^2)*.Cnn(ns[1])
		vd <- ed + (theta^2)
		
		## test statistic, degrees of freedom
		Xsq <- sum(((sa-esa)^2)/vsa) + sum(((div-ed)^2)/vd)
		df <- 1
		
		## compute residuals
		obs <- cbind(sa, div)
		expect <- cbind(esa, ed)
		resids <- (obs-expect)/sqrt(expect)
		
	}
	else {
		stop("This package only supports the 1- or 2-population HKA test.", call. = FALSE)
	}
	
	rez <- list(statistic = Xsq, parameter = df, p.value = 1-pchisq(Xsq, df),
				observed = obs, expected = expect,
				residuals = resids,
				estimates = list(theta = theta, That = Tbig, fhat = fhat))
	return(rez)
	
}

#' Perform the HKA test.
#' 
#' @param x an SFS, or a list thereof (one per locus)
#' @param y another SFS; if \code{x} not a list, this must be provided 
#' @param simulate logical; if \code{TRUE}, get p-values by coalescent simulations instead of analytical approximation
#' @param nsims number of simulations to perform (if \code{simulate == TRUE})
#' @param alpha target significance level for parameter estimates (if \code{simulate == TRUE})
#' @param scale the scaling factor for \eqn{\theta}; 1 = autosome, 3/4 = X chromosome, 1/4 = Y chromosome or mitochondria
#' @param ... ignored
#' 
#' @return a list of class \code{"htest"} with test statistic, degrees of freedom, p-value, and estimates of coalescent parameters 
#' 
#' @references
#' 	Hudson RR, Kreitman M, Aguadé M (1987) A test of neutral molecular evolution based on nucleotide data. Genetics 116: 153-159.
#' 
#' @export
hka_test <- function(x, y = NULL, simulate = FALSE, nsims = 1000, alpha = 0.05, scale = 1.0, ...) {
	
	if (!is.list(x)) {
		xx <- vector("list", 2)
		xx[[1]] <- x
		xx[[2]] <- y
		x <- xx
	}
	x <- lapply(x, noboots)
	x <- lapply(x, matrify)
	if (!same_size(x))
		stop("All input SFS must have same dimensions.", call. = FALSE)
	
	rez <- .hka_test(x, scale)
	if (simulate) {
		nulldist <- .simulate.hka(nsims, dim_sfs(x[[1]]), rez$estimates$theta/scale, rez$estimates$That, rez$estimates$fhat)
		nulldist <- as.data.frame(t(simplify2array(nulldist)))
		rez$null <- nulldist
		rez$p.value <- max(sum(nulldist$statistic > rez$statistic), 1)/nsims
		rez$method <- paste("Tail probabilities obtained by ", nsims, "coalescent simulations")
	}
	else {
		rez$method <- "Tail probabilites from asymptotic chi-squared approximation"
	}
		
	class(rez) <- c(c("htest","hka.test"), class(rez))
	
	return(rez)
	
}

.simulate.hka <- function(nsims, ssz, theta, time, f, path = "ms", ...) {
	
	.harvest.sfs <- function(ltheta) {
		.run.ms(ssz, nsims, ltheta, time, f, path, ...)$sfs
	}
	
	f <- f[1]
	if (is.na(f)) f <- 1
	
	perlocus <- lapply(theta, .harvest.sfs)
	sims <- vector("list", nsims)
	lapply(seq_len(nsims), function(ii) {
		this.sim <- lapply(seq_along(theta), function(f) perlocus[[f]][[ii]])
		#print(this.sim)
		rez <- .hka_test(this.sim)
		params <- do.call("c", rez$estimates)
		c(params, "statistic" = rez$statistic)
	})
	
}

.run.ms <- function(ssz, nsims, theta, time, f, path = "ms", quiet = FALSE, ...) {
	
	n <- sum(ssz)
	cmd <- paste(path, n, nsims, paste("-t", theta))
	if (length(ssz) > 1) {
		cmd <- paste(cmd, paste("-I", length(ssz), paste(ssz, collapse = " ")))
		cmd <- paste(cmd, "-ej", time, 1, 2)
		cmd <- paste(cmd, "-en", time, 1, f)
	}
	
	if (!quiet)
		message(cmd)
	
	ll <- character()
	piper <- pipe(cmd, open = "r")
	while(length(theline <- readLines(piper, 1))) {
		ll <- c(ll, theline)
	}
	sslines <- grep("^segsites", ll)
	close(piper)
	
	ii <- grep("//", ll, fixed = TRUE)
	ss <- as.numeric(gsub("segsites\\: ", "", ll[ii+1]))
	sfs <- vector("list", nsims)
	for (jj in seq_along(ii)) {
		
		if (ss[jj] > 0) {
			start <- ii[jj]+3
			end <- start+sum(ssz)-1
			chroms <- lapply(ll[start:end], function(f) as.integer(strsplit(f, "")[[1]]))
			chroms <- as.data.frame(t(do.call(cbind, chroms)))
			chroms <- .split.ms.pops(chroms, ssz)
			sfs[[jj]] <- .chroms.to.sfs(chroms, ssz)
		}
		else {
			sfs[[jj]] <- matrify(array(0, dim = ssz+1))
		}

	}
	
	return( list(sfs = sfs, segsites = ss) )
	
}

.split.ms.pops <- function(chroms, ssz, ...) {
	
	grouper <- lapply(seq_along(ssz), function(ii) rep(ii, ssz[ii]))
	grouper <- do.call("c", grouper)
	lapply(split(chroms, grouper, drop = FALSE), as.matrix)
	
}

.chroms.to.sfs <- function(chroms, ssz, ...) {

	x <- array(0, dim = ssz+1)
	x <- matrify(x)
	ii <- do.call(cbind, lapply(chroms, colSums))
	for (jj in seq_len(nrow(ii))) {
		idx <- ii[ jj,, drop = FALSE ]+1
		x[ idx ] <- x[ idx ] + 1
	}
	return(x)
	
}

## count number of segregating sites in SFS
.segsites <- function(sfs, ...) {
	
	if (.hasboots(sfs))
		sapply(attr(sfs, "bootstraps"), function(f) sum(mask_corners(f)))
	else
		sum(mask_corners(sfs))
	
}

## helper functions for computing Tajima's D
## see <https://github.com/ANGSD/angsd/blob/master/R/tajima.R>
## n is the number of *chromosomes* in the sample
a1 <- function(n) sum(1/seq_len(n-1))
a2 <- function(n) sum(1/seq_len(n-1)^2)

b1 <- function(n) (n+1)/(3*(n-1))
b2 <- function(n) (2*(n^2+n+3))/(9*n*(n-1))

c1 <- function(n) b1(n)-1/a1(n)
c2 <- function(n) b2(n)-(n+2)/(a1(n)*n)+a2(n)/a1(n)^2

e1 <- function(n) c1(n)/a1(n)
e2 <- function(n) c2(n)/(a1(n)^2+a2(n))

C <- function(n) if (n == 2) 1 else 2*( ((n*a1(n))-2*(n-1))/((n-1)*(n-2)) )
.tajima.denom <- function(n,S) sqrt(e1(n)*S + e2(n)*S*(S-1))

#' Calculate Tajima's D
#' 
#' @param x an SFS
#' @param ... ignored
#' 
#' @references
#' Tajima F (1989) Statistical method for testing the neutral mutation hypothesis by DNA polymorphism. Genetics 105: 437-460.
#' 
#' @export
#' @rdname neutrality_tests
tajimaD <- function(x, ...){
	
	if (ndim_sfs(x) > 1)
		stop("Tajima's D is only defined for a 1d SFS -- try marginalizing first.")
	
	n <- dim_sfs(x)
	pw <- theta_pi(x, persite = FALSE)
	S <- .segsites(x)
	numerator <- pw - S/a1(n)
	return(numerator/.tajima.denom(n,S))
	
}

#' Calculate Fu and Li's D (with outgroup)
#' 
#' @references
#' Fu YX and Li WH (1993) Statistical tests of neutrality of mutations. Genetics 133: 693-709.
#' 
#' @export
#' @rdname neutrality_tests
fuliD <- function(x, ...) {
	
	if (ndim_sfs(x) > 1)
		stop("Fu and Li's D is only defined for a 1d SFS -- try marginalizing first.")
	
	nu <- function(n) 1 + (a1(n)^2)/(a2(n)+(a1(n)^2))*( C(n)-(n+1)/(n-1) )
	uu <- function(n) a1(n) - 1 - nu(n)
	
	n <- dim_sfs(x)
	zeta <- theta_zeta(x, persite = FALSE)
	S <- .segsites(x)
	numerator <- S-zeta*a1(n)
	denominator <- sqrt( S*uu(n) + (S^2)*nu(n) )
	return(numerator/denominator)
	
}

#' Calculate Fu and Li's F (with outgroup)
#' 
#' @references
#' Fu YX and Li WH (1993) Statistical tests of neutrality of mutations. Genetics 133: 693-709.
#' 
#' @export
#' @rdname neutrality_tests
fuliF <- function(x, ...) {
	
	if (ndim_sfs(x) > 1)
		stop("Fu and Li's F is only defined for a 1d SFS -- try marginalizing first.")
	
	vf <- function(n) (C(n) + (2*(n^2+n+3)/(9*n*(n-1))) - 2/(n-1))/(a1(n)^2 + a2(n))
	uf <- function(n) (1 + (n+1)/(3*(n-1)) - ((4*(n+1))/((n-1)^2))*(a1(n+1) - (2*n)/(n+1)))/(a1(n)) - vf(n)
	
	n <- dim_sfs(x)
	pw <- theta_pi(x, persite = FALSE)
	zeta <- theta_zeta(x, persite = FALSE)
	S <- .segsites(x)
	numerator <- pw-zeta
	denominator <- sqrt(S*uf(n)+(S^2)*vf(n))
	return(numerator/denominator)
	
}


#' Calculate Achaz's Y* statistic
#' 
#' @references
#' Achaz G (2008) Testing for neutrality in samples with sequencing errors. Genetics 179: 1409-1424
#' 
#' @export
#' @rdname neutrality_tests
achazY <- function(x, ...) {
	
	if (ndim_sfs(x) > 1)
		stop("Achaz's Y is only defined for a 1d SFS -- try marginalizing first.")
	
	ff <- function(n) (n-2)/(n*(a1(n)-1))
	alpha <- function(n) (ff(n)^2)*(a1(n)-1) + ff(n)*( a1(n)*(4*(n+1)/((n-1)^2)) - 2*(n+1)*(n+2)/(n*(n-1)) ) -
		a1(n)*8*(n-1)/(n*(n-1)^2) + (n^2+n^2+60*n+12)/((3*n^2)*(n-1))
	beta <- function(n) (ff(n)^2)*(a2(n) + a1(n)*(4/((n-1)*(n-2))) - 4/(n-2)) +
		ff(n)*( -a1(n)*(4*(n+2)/(n*(n-1)*(n-2))) - (n^2-3*n^2-16*n+20)/(n*(n-1)*(n-2)) ) +
		a1(n)*( 8/(n*(n-1)*(n-2)) ) + (2*(n^4-n^3-17*n^2-42*n+72))/((9*n^2)*(n-1)*(n-2))
	
	n <- dim_sfs(x)
	y <- mask_singletons(x)
	
	S <- .segsites(y)
	pw <- theta_pi(y, persite = FALSE)
	that <- S/(a1(n)-1)
	thatsq <- S*(S-1)/((a1(n)-1)^2)
	numerator <- pw-ff(n)*S
	denominator <- sqrt(alpha(n)*that + beta(n)*thatsq)
	return(numerator/denominator)
	
}

#' Calculate F_st (genetic differentiation among populations)
#' 
#' @param x an SFS with >1 dimension
#' @param bootstraps logical; if \code{TRUE} (the default), run calculation on any available bootstraps
#' @param ... ignored
#' 
#' @return estimate of \eqn{F_{st}}, via Weir and Cockerham's classic estimator \eqn{\hat{\theta}}
#' 
#' @references
#' Weir BS, Cockerham C (1984) Estimating F-statistics for the analysis of population structure. Evolution 38: 1358-1370.
#' 
#' @export
f_st <- function(x, bootstraps = TRUE, ...) {
	
	.fst <- function(y) {
		## quantities determined only by sample size
		n <- dim_sfs(y)
		r <- ndim_sfs(y)
		nbar <- mean(n)
		nc <- (r*nbar - sum(n^2)/(r*nbar))/(r-1)
		
		## everything else we compute once for bin in the SFS
		idx <- which(!is.na(y), arr.ind = TRUE)
		fst <- numeric(length(y))
		y <- mask_corners(y)
		w <- as.vector(y)/sum(y)
		for (ii in seq_along(y)) {
			cell <- idx[ii,]-1 # absolute allele counts
			pi <- cell/n # convert to frequencies
			pbar <- mean(pi)
			s2 <- var(pi)
			# next two lines may apply to haploids only
			hi <- 1-pi^2
			hbar <- mean(hi)
			# now the fun part
			a <- (nbar/nc)*( s2 - (1/(nbar-1))*(pbar*(1-pbar) - s2*(r-1)/r - hbar/4) )
			b <- (nbar/(nbar-1))*(pbar*(1-pbar) - s2*(r-1)/r - hbar*(2*nbar-1)/(4*nbar))
			c <- hbar/2
			thetahat <- a/(a+b+c)
			fst[ii] <- thetahat
		}
		
		## final value of Fst is that of each bin, times frequency of that bin
		fst <- weighted.mean(fst, w)
		return(fst)
		
	}

	if (ndim_sfs(x) < 2)
		stop("Need an SFS with >= 2 dimensions to calclate F_st.")
		
	if (.hasboots(x) && bootstraps)
		.sapply.boot(x, .fst)
	else
		.fst(x)
	
}