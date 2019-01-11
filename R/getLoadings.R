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
