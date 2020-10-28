#' Principal Component Analysis (PCA) is a very powerful technique that has wide applicability in data science, bioinformatics, and further afield. It was initially developed to analyse large volumes of data in order to tease out the differences/relationships between the logical entities being analysed. It extracts the fundamental structure of the data without the need to build any model to represent it. This 'summary' of the data is arrived at through a process of reduction that can transform the large number of variables into a lesser number that are uncorrelated (i.e. the ‘principal components'), whilst at the same time being capable of easy interpretation on the original data. PCAtools provides functions for data exploration via PCA, and allows the user to generate publication-ready figures. PCA is performed via BiocSingular - users can also identify optimal number of principal component via different metrics, such as elbow method and Horn's parallel analysis, which has relevance for data reduction in single-cell RNA-seq (scRNA-seq) and high dimensional mass cytometry data.
#'
#' @param mat A data-matrix or data-frame containing numerical data only.
#'   Variables are expected to be in the rows and samples in the columns
#'   by default.
#' @param metadata A data-matrix or data-frame containing metadata. This will
#'   be stored in the resulting pca object. Strictly enforced that
#'   rownames(metadata) == colnames(mat).
#' @param center Center the data before performing PCA? Same as prcomp()
#'   'center' parameter.
#' @param scale Scale the data? Same as prcomp() 'scale' parameter.
#' @param rank An integer scalar specifying the number of PCs to retain.
#'   OPTIONAL for an exact SVD, whereby it defaults to all PCs. Otherwise
#'   REQUIRED for approximate SVD methods.
#' @param removeVar Remove this \% of variables based on low variance.
#' @param transposed Is \code{mat} transposed? DEFAULT = FALSE. If set to TRUE,
#'   samples are in the rows and variables are in the columns.
#' @param BSPARAM A \linkS4class{BiocSingularParam} object specifying the 
#'   algorithm to use for the SVD. Defaults to an exact SVD.
#'
#' @details Principal Component Analysis (PCA) is a very powerful technique that has wide applicability in data science, bioinformatics, and further afield. It was initially developed to analyse large volumes of data in order to tease out the differences/relationships between the logical entities being analysed. It extracts the fundamental structure of the data without the need to build any model to represent it. This 'summary' of the data is arrived at through a process of reduction that can transform the large number of variables into a lesser number that are uncorrelated (i.e. the ‘principal components'), whilst at the same time being capable of easy interpretation on the original data. PCAtools provides functions for data exploration via PCA, and allows the user to generate publication-ready figures. PCA is performed via BiocSingular - users can also identify optimal number of principal component via different metrics, such as elbow method and Horn's parallel analysis, which has relevance for data reduction in single-cell RNA-seq (scRNA-seq) and high dimensional mass cytometry data.
#'
#' @return A \code{\link{pca}} object, containing:
#'   \itemize{
#'     \item \code{rotated}, a data frame of the rotated data, i.e., the centred and scaled (
#'       if either or both are requested) input data multiplied by the variable loadings
#'       ('loadings'). This is the same as the 'x' variable returned by prcomp().
#'     \item \code{loadings}, a data frame of variable loadings ('rotation' variable returned
#'       by prcomp()).
#'     \item \code{variance}, a numeric vector of the explained variation for each principal component.
#'     \item \code{sdev}, the standard deviations of the principal components.
#'     \item \code{metadata}, the original metadata
#'     \item \code{xvars}, a character vector of rownames from the input data.
#'     \item \code{yvars}, a character vector of colnames from the input data.
#'     \item \code{components}, a character vector of principal component / eigenvector names.
#'   }
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
#'   getComponents(p)
#'
#'   getVars(p)
#'
#'   getLoadings(p)
#'
#'   screeplot(p)
#'
#'   screeplot(p, hline = 80)
#'
#'   biplot(p)
#'
#'   biplot(p, colby = 'Group', shape = 'Group')
#'
#'   biplot(p, colby = 'Group', colkey = c(A = 'forestgreen', B = 'gold'),
#'     legendPosition = 'right')
#'
#'   biplot(p, colby = 'Group', colkey = c(A='forestgreen', B='gold'),
#'     shape = 'Group', shapekey = c(A=10, B=21), legendPosition = 'bottom')
#'
#'   pairsplot(p, triangle = TRUE)
#'
#'   plotloadings(p, drawConnectors=TRUE)
#'
#'   eigencorplot(p, components = getComponents(p, 1:10),
#'     metavars = c('ESR', 'CRP'))
#'
#' @import BiocParallel
#' @import Matrix
#' @importFrom DelayedArray DelayedArray
#' @importFrom BiocSingular runPCA
#' @importFrom BiocSingular ExactParam
#' @importFrom DelayedMatrixStats colVars
#'
#' @export
pca <- function(
  mat,
  metadata = NULL,
  center = TRUE,
  scale = FALSE,
  rank = NULL, 
  removeVar = NULL,
  transposed = FALSE,
  BSPARAM = ExactParam())
{
  # avoid attempting to coerce S4 matrices into full matrices.
  if (is.data.frame(mat)) {
    mat <- as.matrix(mat)
  }
  if (!transposed) {
    mat <- t(mat)
  }

  # if metadata specified, enforce rule that rownames(metadata) is the
  # same as colnames(mat)
  if (!is.null(metadata)) {
    if(!identical(rownames(mat), rownames(metadata))) {
      stop("'colnames(mat)' is not identical to 'rownames(metadata)'")
    }
  }

  # remove lower portion of variables based on variation
  .center <- if (center) NULL else 0
  vars <- colVars(DelayedArray(mat), center=.center)

  if (!is.null(removeVar)) {
    message('-- removing the lower ', removeVar * 100,
      '% of variables based on variance')
    varorder <- order(vars, decreasing = TRUE)
    keep <- head(varorder, max(1, ncol(mat)*(1-removeVar)))
    mat <- mat[,keep,drop=FALSE]
    vars <- vars[keep]
  }

  # Setting the default rank to all values if Exact.
  if (is.null(rank)) {
    if (is(BSPARAM, "ExactParam")) {
      rank <- min(dim(mat))
    } else {
      stop("'rank' must be specified for approximate PCA methods")
    }
  }

  # perform pca via BiocSingular::runPCA().
  pcaobj <- runPCA(mat,
    center = center,
    scale = scale,
    rank = rank,
    BSPARAM = BSPARAM)

  # Determine the proportion of variance of each component
  # Proportion of variance equals (PC stdev^2) / (total variance)
  if (scale) {
    total.var <- length(vars)
  } else {
    total.var <- sum(vars)
  }
  proportionvar <- (pcaobj$sdev ^ 2) / total.var * 100

  # create a new object
  pcaobj <- list(
    rotated = data.frame(pcaobj$x),
    loadings = data.frame(pcaobj$rotation),
    variance = proportionvar,
    sdev = pcaobj$sdev,
    metadata = metadata,
    xvars = colnames(mat),
    yvars = rownames(mat),
    components = colnames(pcaobj$x)
  )

  rownames(pcaobj$rotated) <- pcaobj$yvars
  rownames(pcaobj$loadings) <- pcaobj$xvars
  names(pcaobj$variance) <- pcaobj$components

  # assign class pca to object
  class(pcaobj) <- 'pca'

  return(pcaobj)
}
