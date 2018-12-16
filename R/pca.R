pca <- function(
  mat,
  metadata = NULL,
  center = TRUE,
  scale = FALSE,
  removeVar = NULL)
{
  mat <- as.data.frame(mat)

  # if metadata specified, enforce rule that rownames(metadata) is the
  # same as colnames(mat)
  if (!is.null(metadata)) {
    if (all((colnames(mat) == rownames(metadata)) == TRUE) == FALSE) {
      stop('Colnames of \'mat\' object must equal and be in the same',
        ' order as the rownames of metadata')
    }
  }

  # remove lower portion of variables based on variation
  if (!is.null(removeVar)) {
    message('-- removing the lower ', removeVar * 100,
      '% of variables based on variance')
    vars <- apply(mat, 1, function(x) var(x))
    varorder <- order(vars, decreasing = FALSE)
    exclude <- varorder[seq_len(nrow(mat)*removeVar)]
    mat <- mat[-exclude,]
  }

  # perform pca via prcomp()
  pcaobj <- prcomp(t(mat),
    center = center,
    scale. = scale,
    retx = TRUE,
    tol = NULL,
    rank. = NULL)

  # Determine the proportion of variance of each component
  # Proportion of variance equals (PC stdev^2) / (sum all PCs stdev^2)
  proportionvar <- ((pcaobj$sdev ^ 2) / (sum(pcaobj$sdev ^ 2))) * 100

  # create a new object
  pcaobj <- list(
    rotated = data.frame(pcaobj$x),
    loadings = data.frame(pcaobj$rotation),
    variance = proportionvar,
    metadata = metadata,
    xvars = rownames(mat),
    yvars = colnames(mat),
    components = colnames(pcaobj$x)
  )

  names(pcaobj$variance) <- pcaobj$components

  # assign class pca to object
  class(pcaobj) <- 'pca'

  return(pcaobj)
}
