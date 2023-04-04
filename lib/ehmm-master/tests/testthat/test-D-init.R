context("D: Initialization algorithm works")
source("utils.R")
#kullback-leibler divergence of two negative multinomials with marginal means mus1 and mus2
#r is the same for both distributions
#that's a non-symmetric distance, it's the distance from mus1 to mus2
#vectorized on the second argument
kldiv <- function(mus1, mus2, r=Inf){
	if (!is.matrix(mus2)) dim(mus2) <- c(length(mus2), 1)
	ratioTerm <- mus1*log(mus1/mus2)
	ratioTerm[mus1==0] <- 0
	ratioTerm[mus1>0&mus2==0] <- Inf
	if (r==Inf){
		return(colSums(mus2-mus1+ratioTerm))
	} else {
		mu1 <- sum(mus1)
		mu2 <- colSums(mus2)
		return((mu1+r)*log((mu2+r)/(mu1+r)) + colSums(ratioTerm))
	}
}

kldmat <- function(mus, r=Inf){
	0.5*(t(apply(mus, 2, kldiv, mus2=mus, r=r)) + apply(mus, 2, kldiv, mus2=mus, r=r))
}

test_that("KL divergence works",{
	#each column represents a neg multinom
	mus <- matrix(
		c(0,0,
			0,3,
			1,12,
			3,3), nrow=2)
	
	#test with a finite r
	r <- 2.5
	dmat <- KL_dist_mat(mus, r)
	expect_equal(diag(dmat), rep(0, ncol(dmat)))
	expect_equal(dmat, t(dmat))
	expect_true(all(dmat>=0))
	expect_true(all(dmat[c(1,2), c(3,4)]==Inf))
	
	#compare with slow algo
	their_dmat <- kldmat(mus, r)
	expect_equal(their_dmat, dmat)
	
	
	#test more random pairs
	mus <- matrix(abs(rnorm(10*100)), nrow=10)
	
	#test with r=Inf
	r <- Inf
	dmat <- KL_dist_mat(mus, r)
	expect_equal(diag(dmat), rep(0, ncol(dmat)))
	expect_equal(dmat, t(dmat))
	expect_true(all(dmat>=0))
	
	#compare with slow algo
	their_dmat <- kldmat(mus, r)
	expect_equal(their_dmat, dmat)
	
	#test with a finite r
	r <- 2.5
	dmat <- KL_dist_mat(mus, r)
	expect_equal(diag(dmat), rep(0, ncol(dmat)))
	expect_equal(dmat, t(dmat))
	expect_true(all(dmat>=0))
	
	#compare with slow algo
	their_dmat <- kldmat(mus, r)
	expect_equal(their_dmat, dmat)
	
	#test with many threads
	expect_equal(KL_dist_mat(mus, Inf), KL_dist_mat(mus, Inf, nthreads=10))
	expect_equal(KL_dist_mat(mus, 2.5), KL_dist_mat(mus, 2.5, nthreads=10))
	
	
})

test_that("tpca works",{
	mat <- matrix(rnorm(10*1000), nrow=10)
	their_cov <- (mat %*% t(mat))/(ncol(mat)-1)
	my_cov <- rowdotprod(mat, besselCorr=T)
	expect_equal(their_cov, my_cov)
	
	my_cov_mc <- rowdotprod(mat, besselCorr=T, nthreads=10)
	expect_equal(their_cov, my_cov_mc)
	
	their_pca <- prcomp(t(mat))
	#adjust 'rotation' and 'x' so that they use the same conventions
	#sign of an eigenvector (convention): sign of the element with highest absolute value
	vsign <- apply(their_pca$rotation, 2, function(col) sign(col[which.max(abs(col))]))
	#make all eigenvectors with positive sign
	their_pca$rotation <- their_pca$rotation*vsign[col(their_pca$rotation)]
	#adjust the x matrix accordingly
	their_pca$x <- their_pca$x*vsign[col(their_pca$x)]
	#one thread
	my_pca <- tpca(mat)
	expect_equal(their_pca$sdev, my_pca$sdev)
	expect_equal(their_pca$rotation, my_pca$rotation)
	expect_equal(their_pca$center, my_pca$center)
	expect_equal(t(their_pca$x), my_pca$tx)
	#many threads
	my_pca_mc <- tpca(mat, nthreads=10)
	expect_equal(their_pca$sdev, my_pca_mc$sdev)
	expect_equal(their_pca$rotation, my_pca_mc$rotation)
	expect_equal(their_pca$center, my_pca_mc$center)
	expect_equal(t(their_pca$x), my_pca_mc$tx)
})

#replace v with a rank value, such that all
#possible ranks are from 0 to nsplit-1 and
#each rank has a similar number of members
their_splitAxes <- function(v, nsplit){
	o <- order(v)
	io <- integer(length(v))
	io[o] <- 1:length(v)
	mult <- nsplit/length(v)
	floor(mult*(io-1))
}

their_clustAvg <- function(counts, coords, clusts){
	clustPerAxis <- max(sapply(clusts, max)) + 1
	clustAvg1Axis <- function(counts, coord, clust){
		coord2 <- clust[coord+1]
		mus <- sapply(1:clustPerAxis, function(i){
			rowMeans(counts[, coord2==(i-1)])
		})
		sizes <- table(factor(coord2, levels=(0:(clustPerAxis-1))))
		mus[,sizes==0] <- 0
		list(mus=mus, sizes=sizes)
	}
	clustList <- lapply(1:length(clusts), function(i){
		clustAvg1Axis(counts, coords[i,], clusts[[i]])
	})
	
	musList <- lapply(clustList, function(clust) clust$mus)
	sizesList <- lapply(clustList, function(clust) clust$sizes)
	mus <- do.call(cbind, musList)
	sizes <- unlist(sizesList); names(sizes) <- NULL
	
	list(mus=mus, sizes=sizes)
}

their_fillPost <- function(coords, clusts, nclust){
	post <- matrix(0, nclust, ncol(coords))
	for (rw in 1:nrow(coords)){
		rclusts <- clusts[[rw]]
		for (cl in 1:ncol(coords)){
			clust <- rclusts[coords[rw,cl]+1]+1
			post[clust,cl] <- post[clust,cl]+1
		}
	}
	post
}

test_that("splitAxes, clustAvg and fillPosteriors",{
	nsplit <- 17
	nclust <- 10
	mat <- matrix(rnorm(10*1000), nrow=10)
	
	#splitAxes
	their_splits <- t(apply(mat, 1, their_splitAxes, nsplit=nsplit))
	my_splits <- splitAxes(mat, nsplit)
	expect_equal(their_splits, my_splits)
	
	my_splits_mc <- splitAxes(mat, nsplit, nthreads=10)
	expect_equal(their_splits, my_splits_mc)
	
	imat <- matrix(rnbinom(10*1000, mu=4, size=1), nrow=10)
	their_splits_int <- t(apply(imat, 1, their_splitAxes, nsplit=nsplit))
	my_splits_int <- splitAxesInt(imat, nsplit)
	expect_equal(their_splits_int, my_splits_int)
	
	my_splits_int_mc <- splitAxesInt(imat, nsplit, nthreads=10)
	expect_equal(their_splits_int, my_splits_int_mc)
	
	#clustAvg
	
	#creating some fake labels
	labs <- lapply(1:nrow(mat), function(i){
		sample(nclust, nsplit, replace=T)-1
	})
	
	their_clusts <- their_clustAvg(imat, their_splits_int, labs)
	my_clusts <- clusterAverages2(imat, their_splits_int, labs)
	expect_equal(their_clusts, my_clusts)
	
	my_clusts_mc <- clusterAverages2(imat, their_splits_int, labs, nthreads=10)
	expect_equal(their_clusts, my_clusts_mc)
	
	#fillPosteriors
	their_post <- their_fillPost(their_splits_int, labs, nclust)
	my_post <- fillPosteriors(their_splits_int, labs, nclust)
	expect_equal(their_post, my_post)
	
	my_post_mc <- fillPosteriors(their_splits_int, labs, nclust, nthreads=10)
	expect_equal(their_post, my_post_mc)
	
})

test_that("the whole thing runs",{
	mat <- matrix(rnbinom(10*1000, mu=4, size=1), nrow=10)
	nclust <- 5
	nthreads <- 1
	nbtype <- "dep"
	nlevs <- c(10,10,100,100)
	axes <- c("counts", "pca", "counts", "pca")
	
	inits <- lapply(1:4, function(i){
		initAlgo(mat, nclust, nlev=nlevs[i], axes=axes[i], nthreads=nthreads, nbtype=nbtype)
	})
	
	for (i in seq_along(inits)){
		init <- inits[[i]]
		nlev <- nlevs[i]
		expect_true(modelsAreOk(init$models, nclust, nrow(mat), nbtype))
		expect_equal(length(init$mix_coeff), nclust)
		expect_true(all(is.finite(init$mix_coeff)))
	}
	
	nthreads <- 10
	inits_mc <- lapply(1:4, function(i){
		initAlgo(mat, nclust, nlev=nlevs[i], axes=axes[i], nthreads=nthreads, nbtype=nbtype)
	})
	
	for (i in seq_along(inits)){
		expect_equal(inits[[i]], inits_mc[[i]], tol=1e-6)
	}
})


test_that("getUniqueSeeds works",{
	mat <- matrix(c(0,0,0,1,0,2,0,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0), nrow=2)
	s1 <- getUniqueSeeds(mat, 1);
	expect_true(s1 > 0 && s1 <= ncol(mat))
	expect_error(getUniqueSeeds(mat, 5))
	s4 <- getUniqueSeeds(mat, 4)
	expect_equal(length(s4), 4)
	expect_true(!any(duplicated(s4)))
	expect_true(all(c(2,3,4) %in% s4))
	
	rmat <- matrix(nrow=4, rpois(4*10000, lambda=0.01))
	rmat <- cbind(rmat, c(0,0,0,1), c(0,0,1,0), c(0,0,1,1), c(0,1,0,0))
	
	s4 <- getUniqueSeeds(rmat, 4)
	expect_equal(length(s4), 4)
	expect_true(!any(duplicated(s4)))
	expect_true(!any(duplicated(rmat[,s4])))
	
})

