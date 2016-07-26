# sfsr
An `R` package for analysing population-genetic data represented as the site frequency spectrum (SFS)

## Background
The *site frequency spectrum* (SFS; also known as the allele frequency spectrum) is the histogram over allele counts in a population.  More formally, for a population consisting of $N$ chromosomes, the SFS is the vector $\mathbf{\zeta} = \zeta_0, \zeta_1, \dots, \zeta_N$ whose entries correspond to the number of sites at which $0, 1, \dots, N$ copies of the derived allele are segregating in the population.  When the polarity (ancestral vs derived) of alleles is unknown, the SFS is said to be "folded;" the folded SFS is often denoted instead $\mathbf{\nu} = \nu_0, \dots, \nu_m$ (with $m = \ceil{\frac{N}{2}}$) and its entries are the count of the minor allele.

The SFS can be generalized to higher dimensions when more than one population is considered; in this case we may call it the joint SFS.  In this case each entry $\mathbf{\zeta}_{i,j,\dots}$ gives the count of alleles with frequency $i$ in population 1, $j$ in population 2, etc.  Given a joint SFS, the single-population SFS can always be obtained by "marginalizing" -- just summing over all other dimensions.

The SFS provides a rich summary of polymorphism data.  Departures in the shape of the SFS are informative for demography and/or selection.  Most common summary statistics for polymorphism data both within and between populations can be computed directly from the SFS.

## Rationale
Although the SFS can be easily tabulated from discrete genotypes, it may be advantageous to treat it as a parameter to be estimated from data under some model of uncertainty and technical error.  Such methods are implemented in software like `ANGSD`, but I found few resources for manipulating and analyzing the resulting SFS.  Hence the `sfsr` package.

At the moment the package supports only unfolded SFS -- those for which alleles can be polarized without ambiguity as ancestral vs derived.  This happens to be the case for problems I work on, and gets around the data-dependency of the "minor allele" designation.

## Features

* File I/O
	* read `ANGSD`-style SFS from disk (`read_sfs()`)
* Manipulation and transformation of SFS
	* marginalization (`margnialize_sfs()`)
	* changes of ancestral-vs-derived polarity (`repolarize_sfs()`)
* Common within-population summary statistics
	* Estimators of $\theta$
		* Tajima's pairwise, aka $\pi$ (`theta_pi()`)
		* Watterson's (`theta_w()`)
		* from singleton count (`theta_zeta()`)
	* Neutrality statistics
		* Tajima's $D$ (`tajimaD()`)
		* Fu and LI's $D$ (`fuliD()`)
		* Fu and LI's $F$ (`fuliF()`)
		* Achaz's $Y$, robust to sequencing errors (`achazY()`)
* Between-population summary statistics
	* fixed vs polymorhic sites (`fixpoly()`)
	* average divergence (`d_xy()`)
	* Weir and Cockerham's $F_{st}$ (`f_st()`)
* Hudson-Kreitman-Aguade test (`hka_test()`)

Because bootstrap resampling from the SFS is a straightforward way to estimate uncertainty, SFS as reprsented by `sfsr` may carry bootstrap replicates as an attribute.  In this case most summary statistics are automatically calculated across the bootstraps for convenient estimation of standard errors.
