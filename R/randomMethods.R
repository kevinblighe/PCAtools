#' Choosing PCs with the Marchenko-Pastur limit
#'
#' Use the Marchenko-Pastur limit to choose the number of top PCs to retain.
#'
#' @param x The data matrix used for the PCA, containing variables in rows and observations in columns.
#' Ignored if \code{dim} is supplied.
#' @param .dim An integer vector containing the dimensions of the data matrix used for PCA.
#' The first element should contain the number of variables and the second element should contain the number of observations.
#' @param var.explained A numeric vector containing the variance explained by successive PCs.
#' This should be sorted in decreasing order.
#' Note that this should be the variance explained, NOT the percentage of variance explained!
#' @param noise Numeric scalar specifying the variance of the random noise.
#'
#' @details
#' For a random matrix with i.i.d. values, the Marchenko-Pastur (MP) limit defines the maximum eigenvalue.
#' Let us assume that \code{x} is the sum of some low-rank truth and some i.i.d. random matrix with variance \code{noise}.
#' We can use the MP limit to determine the maximum variance that could be explained by a fully random PC;
#' all PCs that explain more variance are thus likely to contain real structure and should be retained.
#'
#' Of course, this has some obvious caveats such as the unrealistic i.i.d. assumption and the need to estimate \code{noise}.
#' Moreover, PCs below the MP limit are not necessarily uninformative or lacking structure;
#' it is just that their variance explained does not match the most extreme case that random noise has to offer.
#'
#' @return
#' An integer scalar specifying the number of PCs with variance explained beyond the MP limit.
#' The limit itself is returned in the attributes.
#'
#' @author Aaron Lun
#'
#' @seealso
#' \code{\link{chooseGavishDonoho}}, \code{\link{parallelPCA}} and \code{\link{findElbowPoint}}, 
#' for other approaches to choosing the number of PCs.
#'
#' @examples
#' truth <- matrix(rnorm(1000), nrow=100)
#' truth <- truth[,sample(ncol(truth), 1000, replace=TRUE)]
#' obs <- truth + rnorm(length(truth), sd=2)
#' 
#' # Note, we need the variance explained, NOT the percentage
#' # of variance explained! 
#' pcs <- pca(obs)
#' chooseMarchenkoPastur(obs, var.explained=pcs$sdev^2, noise=4)
#' 
#' @export
chooseMarchenkoPastur <- function(x, .dim=dim(x), var.explained, noise) {
    # Using the Marchenko-Pastur limit on the _eigenvalues_ of the noise matrix with variance 'noise'.
    # (See https://www.wolfram.com/language/11/random-matrices/marchenko-pastur-distribution.html?product=mathematica)
    limit <- RMTstat::qmp(1, ndf=.dim[2] - 1, pdim=.dim[1], var=noise)
    marchenko <- sum(var.explained > limit)

    output <- max(1L, marchenko)
    attr(output, "limit") <- limit
    output
}

#' Choosing PCs with the Gavish-Donoho method
#'
#' Use the Gavish-Donoho method to determine the optimal number of PCs to retain.
#'
#' @inheritParams chooseMarchenkoPastur
#'
#' @details
#' Assuming that \code{x} is the sum of some low-rank truth and some i.i.d. random matrix with variance \code{noise},
#' the Gavish-Donoho method defines a threshold on the singular values that minimizes the reconstruction error from the PCs.
#' This provides a mathematical definition of the \dQuote{optimal} choice of the number of PCs for a given matrix,
#' though it depends on both the i.i.d. assumption and an estimate for \code{noise}.
#'
#' @return
#' An integer scalar specifying the number of PCs to retain.
#' The effective limit on the variance explained is returned in the attributes.
#'
#' @author Aaron Lun
#'
#' @examples
#' truth <- matrix(rnorm(1000), nrow=100)
#' truth <- truth[,sample(ncol(truth), 1000, replace=TRUE)]
#' obs <- truth + rnorm(length(truth), sd=2)
#' 
#' # Note, we need the variance explained, NOT the percentage
#' # of variance explained! 
#' pcs <- pca(obs)
#' chooseGavishDonoho(obs, var.explained=pcs$sdev^2, noise=4)
#'
#' @seealso
#' \code{\link{chooseMarchenkoPastur}}, \code{\link{parallelPCA}} and \code{\link{findElbowPoint}}, 
#' for other approaches to choosing the number of PCs.
#' 
#' @export
chooseGavishDonoho <- function(x, .dim=dim(x), var.explained, noise) {
    m <- min(.dim)
    n <- max(.dim)
    beta <- m/n

    # Equation 11 of the Gavish-Donoho paper.
    lambda <- sqrt( 2 * (beta + 1) + (8 * beta) / ( beta + 1 + sqrt(beta^2 + 14 * beta + 1) ) )

    # Equation 3, slightly reorganized to have the variance explained on the LHS
    # instead of the singular values 'D', where D = [ var.explained * (.dim[2] - 1) ].
    limit <- lambda^2 * noise * n/(.dim[2] - 1)
    gv <- sum(var.explained > limit)

    output <- max(1L, gv)
    attr(output, "limit") <- limit
    output
}
