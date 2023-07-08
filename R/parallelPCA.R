#' Perform Horn's parallel analysis to choose the number of principal components to retain.
#'
#' @param mat A numeric matrix where rows correspond to variables and columns correspond to samples.
#' @param max.rank Integer scalar specifying the maximum number of PCs to retain.
#' @param ... Further arguments to pass to \code{\link{pca}}.
#' @param niters Integer scalar specifying the number of iterations to use for the parallel analysis.
#' @param threshold Numeric scalar representing the \dQuote{p-value} threshold above which PCs are to be ignored.
#' @param transposed Logical scalar indicating whether \code{mat} is transposed, i.e., rows are samples and columns are variables.
#' @param BSPARAM A \linkS4class{BiocSingularParam} object specifying the algorithm to use for PCA.
#' @param BPPARAM A \linkS4class{BiocParallelParam} object specifying how the iterations should be paralellized.
#'
#' @details Horn's parallel analysis involves shuffling observations within each row of
#'   \code{x} to create a permuted matrix.  PCA is performed on the permuted matrix
#'   to obtain the percentage of variance explained under a random null hypothesis.
#'   This is repeated over several iterations to obtain a distribution of curves on
#'   the scree plot.
#'
#'   For each PC, the \dQuote{p-value} (for want of a better word) is defined as the
#'   proportion of iterations where the variance explained at that PC is greater
#'   than that observed with the original matrix. The number of PCs to retain is
#'   defined as the last PC where the p-value is below \code{threshold}. This aims
#'   to retain all PCs that explain \dQuote{significantly} more variance than
#'   expected by chance.
#'
#'   This function can be sped up by specifying \code{BSPARAM=IrlbaParam()} or
#'   similar, to use approximate strategies for performing the PCA.  Another option
#'   is to set \code{BPPARAM} to perform the iterations in parallel.
#'
#' @return A list is returned, containing:
#'   \itemize{
#'     \item \code{original}, the output from running \code{\link{pca}} on \code{mat} with the specified arguments.
#'     \item \code{permuted}, a matrix of variance explained from randomly permuted matrices. 
#'       Each column corresponds to a single permutated matrix, while each row corresponds to successive principal components.
#'     \item  \code{n}, the estimated number of principal components to retain.
#'   }
#'
#' @author Aaron Lun
#'
#' @examples
#'   # Mocking up some data.
#'   ngenes <- 1000
#'   means <- 2^runif(ngenes, 6, 10)
#'   dispersions <- 10/means + 0.2
#'   nsamples <- 50
#'   counts <- matrix(rnbinom(ngenes*nsamples, mu=means, 
#'     size=1/dispersions), ncol=nsamples)
#'
#'   # Choosing the number of PCs
#'   lcounts <- log2(counts + 1)
#'   output <- parallelPCA(lcounts)
#'   output$n
#'
#' @importFrom dqrng generateSeedVectors
#' @importFrom BiocParallel bpmapply
#' @importFrom BiocParallel SerialParam
#'
#' @export
parallelPCA <- function(mat, max.rank=100, ..., niters=50, threshold=0.1, 
  transposed=FALSE, BSPARAM=ExactParam(), BPPARAM=SerialParam()) 
{
  if (!transposed) {
    mat <- t(mat)
  }
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

#' @importFrom beachmat initializeCpp
.parallel_PA <- function(mat, max.rank, ..., seed, stream, BSPARAM)
# Function for use in bplapply, defined here to automatically 
# take advantage of the namespace when using snowParam. Note 
# that shuffle_matrix needs samples in columns, hence the t()
# above and transposed=TRUE in the pca() calls.
{
  ptr <- initializeCpp(mat)
  re.y <- shuffle_matrix(ptr, seed, stream)
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
