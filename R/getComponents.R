getComponents <- function(
  pcaobj,
  components = NULL)
{
  if (!is.null(components)) {
    return(pcaobj$components[components])
  } else if (is.null(components)) {
    return(pcaobj$components)
  }
}
