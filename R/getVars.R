getVars <- function(
  pcaobj,
  components = NULL)
{
  if (!is.null(components)) {
    return(pcaobj$variance[components])
  } else if (is.null(components)) {
    return(pcaobj$variance)
  }
}
