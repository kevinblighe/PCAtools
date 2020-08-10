#' Return component loadings for principal components from an object of class 'pca'.
#'
#' @param pcaobj Object of class 'pca' created by pca().
#' @param components Indices of the principal components whose component loadings
#'   will be returned. If NULL, all PC names will be returned.
#'
#' @details Return component loadings for principal components from an object of class 'pca'.
#'
#' @return A \code{\link{data.frame}} object.
#'
#' @author Kevin Blighe <kevin@clinicalbioinformatics.co.uk>
#'
#' @examples
#'   options(scipen=10)
#'   options(digits=6)
#'
#'   col <- 20
#'   row <- 20000
#'   mat1 <- matrix(
#'     rexp(col*row, rate = 0.1),
#'     ncol = col)
#'   rownames(mat1) <- paste0('gene', 1:nrow(mat1))
#'   colnames(mat1) <- paste0('sample', 1:ncol(mat1))
#'
#'   mat2 <- matrix(
#'     rexp(col*row, rate = 0.1),
#'     ncol = col)
#'   rownames(mat2) <- paste0('gene', 1:nrow(mat2))
#'   colnames(mat2) <- paste0('sample', (ncol(mat1)+1):(ncol(mat1)+ncol(mat2)))
#'
#'   mat <- cbind(mat1, mat2)
#'
#'   metadata <- data.frame(row.names = colnames(mat))
#'   metadata$Group <- rep(NA, ncol(mat))
#'   metadata$Group[seq(1,40,2)] <- 'A'
#'   metadata$Group[seq(2,40,2)] <- 'B'
#'   metadata$CRP <- sample.int(100, size=ncol(mat), replace=TRUE)
#'   metadata$ESR <- sample.int(100, size=ncol(mat), replace=TRUE)
#'
#'   p <- pca(mat, metadata = metadata, removeVar = 0.1)
#'
#'   getLoadings(p)
#' 
#' @export
getLoadings <- function(
  pcaobj,
  components = NULL)
{
  if (!is.null(components)) {
    return(pcaobj$loadings[components])
  } else if (is.null(components)) {
    return(pcaobj$loadings)
  }
}
