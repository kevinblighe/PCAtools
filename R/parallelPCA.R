parallelPCA <- function(mat, max.rank=100, ..., niters=50, threshold=0.1, 
  BSPARAM=ExactParam(), BPPARAM=SerialParam()) 
{
  mat <- t(mat)
  original <- pca(mat, rank=max.rank, ..., transposed=TRUE, BSPARAM=BSPARAM)
  original.s2 <- original$variance

  # Running across permutations.
  pcg.states <- .setup_pcg_state(niters)
  permuted <- bpmapply(FUN=.parallel_PA,
    seed=pcg.states$seeds[[1]], stream=pcg.states$streams[[1]],
    MoreArgs=list(mat=mat, ..., max.rank=max.rank, BSPARAM=BSPARAM), 
    BPPARAM=BPPARAM, SIMPLIFY=FALSE, USE.NAMES=FALSE)
  permutations <- do.call(cbind, permuted)

  # Figuring out where the original drops to "within range" of permuted.
  prop <- rowMeans(permutations >= original.s2)
  above <- prop > threshold
  if (!any(above)) {
    npcs <- length(above) 
  } else {
    npcs <- min(which(above)) - 1L
  }

  list(original=original, permuted=permutations, n=npcs)
}

.parallel_PA <- function(mat, max.rank, ..., seed, stream, BSPARAM)
# Function for use in bplapply, defined here to automatically 
# take advantage of the namespace when using snowParam. Note 
# that shuffle_matrix needs samples in columns, hence the t()
# above and transposed=TRUE in the pca() calls.
{
  re.y <- shuffle_matrix(mat, seed, stream)
  out <- pca(re.y, rank=max.rank, ..., transposed=TRUE, BSPARAM=BSPARAM)
  out$variance
}

.setup_pcg_state <- function(per.core) 
# Sets up seeds in the serial component to ensure that the same 
# stream of random numbers is used in each worker, regardless 
# of the parallelization scheme.
{
  seeds <- streams <- vector("list", length(per.core))
  last <- 0L
  for (i in seq_along(per.core)) {
    N <- per.core[i]
    seeds[[i]] <- generateSeedVectors(N, nwords=2)
    streams[[i]] <- last + seq_len(N)
    last <- last + N
  }
  list(seeds=seeds, streams=streams)
}
