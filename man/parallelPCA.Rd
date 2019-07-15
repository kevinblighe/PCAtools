\name{parallelPCA}
\alias{parallelPCA}

\title{Horn's parallel analysis}
\description{Perform Horn's parallel analysis to choose the number of principal components to retain.}

\usage{
parallelPCA(mat, max.rank=100, ..., niters=50, threshold=0.1, 
  BSPARAM=EmatactParam(), BPPARAM=SerialParam())
}

\arguments{
\item{mat}{A numeric matrix where rows correspond to variables and columns correspond to samples.}
\item{max.rank}{Integer scalar specifying the maximum number of PCs to retain.}
\item{...}{Further arguments to pass to \code{\link{pca}}.}
\item{niters}{Integer scalar specifying the number of iterations to use for the parallel analysis.}
\item{threshold}{Numeric scalar representing the \dQuote{p-value} threshold above which PCs are to be ignored.}
\item{BSPARAM}{A \linkS4class{BiocSingularParam} object specifying the algorithm to use for PCA.}
\item{BPPARAM}{A \linkS4class{BiocParallelParam} object specifying how the iterations should be paralellized.}
}

\value{
A list is returned, containing:
\itemize{
  \item \code{original}, the output from running \code{\link{pca}} on \code{mat} with the specified arguments.
  \item \code{permuted}, a matrix of variance explained from randomly permuted matrices. 
  Each column corresponds to a single permutated matrix, while each row corresponds to successive principal components.
 \item  \code{n}, the estimated number of principal components to retain.
}
}

\author{
Aaron Lun
}

\details{
Horn's parallel analysis involves shuffling observations within each row of
\code{x} to create a permuted matrix.  PCA is performed on the permuted matrix
to obtain the percentage of variance explained under a random null hypothesis.
This is repeated over several iterations to obtain a distribution of curves on
the scree plot.

For each PC, the \dQuote{p-value} (for want of a better word) is defined as the
proportion of iterations where the variance explained at that PC is greater
than that observed with the original matrix. The number of PCs to retain is
defined as the last PC where the p-value is below \code{threshold}. This aims
to retain all PCs that explain \dQuote{significantly} more variance than
expected by chance.

This function can be sped up by specifying \code{BSPARAM=IrlbaParam()} or
similar, to use approximate strategies for performing the PCA.  Another option
is to set \code{BPPARAM} to perform the iterations in parallel.
}

\references{
Horn JL. (1965). 
A rationale and test for the number of factors in factor analysis.
\emph{Psychometrika} 30(2), 179-185. 

Buja A and Eyuboglu N (1992).
Remarks on Parallel Analysis.
\emph{Multivariate Behav. Res.} 27, 509-40.
}

\examples{
# Mocking up some data.
ngenes <- 1000
means <- 2^runif(ngenes, 6, 10)
dispersions <- 10/means + 0.2
nsamples <- 50
counts <- matrix(rnbinom(ngenes*nsamples, mu=means, 
    size=1/dispersions), ncol=nsamples)

# Choosing the number of PCs
lcounts <- log2(counts + 1)
output <- parallelPCA(lcounts)
output$n
}
