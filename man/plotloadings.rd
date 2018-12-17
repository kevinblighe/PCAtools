\name{plotloadings}

\alias{plotloadings}

\title{PCAtools: everything Principal Components Analysis}

\description{Plot the component loadings for selected principal components / eigenvectors and label variables driving variation along these.}

\usage{
  plotloadings(pcaobj,
  components = getComponents(pcaobj, seq(1, 5)),
  rangeRetain = 0.05,
  absolute = FALSE,
  col = c('gold', 'white', 'royalblue'),
  colMidpoint = 0,
  shape = 21,
  shapeSizeRange = c(10, 10),
  legendPosition = 'top',
  legendLabSize = 10,
  legendIconSize = 3.0,
  xlim = NULL,
  ylim = NULL,
  labSize = 2.0,
  labhjust = 1.5,
  labvjust = 0,
  drawConnectors = TRUE,
  positionConnectors = 'right',
  widthConnectors = 0.5,
  typeConnectors = 'closed',
  endsConnectors = 'first',
  lengthConnectors = unit(0.01, 'npc'),
  colConnectors = 'grey50',
  xlab = 'Principal component',
  xlabAngle = 0,
  xlabhjust = 0.5,
  xlabvjust = 0.5,
  ylab = 'Component loading',
  ylabAngle = 0,
  ylabhjust = 0.5,
  ylabvjust = 0.5,
  axisLabSize = 16,
  title = '',
  subtitle = '',
  caption = '',
  titleLabSize = 16,
  subtitleLabSize = 12,
  captionLabSize = 12,
  hline = c(0),
  hlineType = 'longdash',
  hlineCol = 'black',
  hlineWidth = 0.4,
  vline = NULL,
  vlineType = 'longdash',
  vlineCol = 'black',
  vlineWidth = 0.4,
  gridlines.major = TRUE,
  gridlines.minor = TRUE,
  borderWidth = 0.8,
  borderColour = 'black',
  returnPlot = TRUE)
}

\arguments{
  \item{pcaobj}{Object of class 'pca' created by pca(). REQUIRED.}
  \item{components}{The principal components to be included in the plot.
  DEFAULT = getComponents(pcaobj, seq(1, 5)). OPTIONAL.}
  \item{rangeRetain}{Cut-off value for retaining variables. The function
  will look across each specified principal component and retain the variables
  that fall within this top/bottom fraction of the loadings range. DEFAULT
  = 0.05. OPTIONAL.}
  \item{absolute}{Plot absolute loadings? (TRUE / FALSE). DEFAULT = FALSE.
  OPTIONAL.}
  \item{col}{Colours used for generation of fill gradient according to
  loadings values. Can be 2 or 3 colours. DEFAULT =
  c('gold', 'white', 'royalblue'). OPTIONAL.}
  \item{colMidpoint}{Mid-point (loading) for the colour range. DEFAULT = 0.
  OPTIONAL.}
  \item{shape}{Shape of the plotted points. DEFAULT = 21. OPTIONAL.}
  \item{shapeSizeRange}{Size range for the plotted points (min, max). DEFAULT
  = c(10, 10). OPTIONAL.}
  \item{legendPosition}{Position of legend ('top', 'bottom', 'left', 'right',
  'none'). DEFAULT = 'top'. OPTIONAL.}
  \item{legendLabSize}{Size of plot legend text. DEFAULT = 10. OPTIONAL.}
  \item{legendIconSize}{Size of plot legend icons / symbols. DEFAULT = 3.0.
  OPTIONAL.}
  \item{xlim}{Limits of the x-axis. DEFAULT = NULL. OPTIONAL.}
  \item{ylim}{Limits of the y-axis. DEFAULT = NULL. OPTIONAL.}
  \item{labSize}{Size of labels. DEFAULT = 2.0. OPTIONAL.}
  \item{labhjust}{Horizontal adjustment of label. DEFAULT = 1.5. OPTIONAL.}
  \item{labvjust}{Vertical adjustment of label. DEFAULT = 0. OPTIONAL.}
  \item{drawConnectors}{Fit labels onto plot and connect to their respective
  points by line connectors (TRUE/FALSE). DEFAULT = FALSE. OPTIONAL.}
  \item{positionConnectors}{Position of the connectors and their labels with
  respect to the plotted points ('left', 'right'). DEFAULT = 'right'.
  OPTIONAL.}
  \item{widthConnectors}{Line width of connectors. DEFAULT = 0.5. OPTIONAL.}
  \item{typeConnectors}{Have the arrow head open or filled ('closed')?
  ('open', 'closed'). DEFAULT = 'closed'. OPTIONAL.}
  \item{endsConnectors}{Which end of connectors to draw arrow head? ('last',
  'first', 'both'). DEFAULT = 'first'. OPTIONAL.}
  \item{lengthConnectors}{Length of the connectors. DEFAULT =
  unit(0.01, 'npc'). OPTIONAL}
  \item{colConnectors}{Line colour of connectors. DEFAULT = 'grey50'. OPTIONAL.}
  \item{xlab}{Label for x-axis. DEFAULT = 'Principal component'. OPTIONAL.}
  \item{xlabAngle}{Rotation angle of x-axis labels. DEFAULT = 0. OPTIONAL.}
  \item{xlabhjust}{Horizontal adjustment of x-axis labels. DEFAULT = 0.5.
  OPTIONAL.}
  \item{xlabvjust}{Vertical adjustment of x-axis labels. DEFAULT = 0.5.
  OPTIONAL.}
  \item{ylab}{Label for y-axis. DEFAULT = 'Component loading'. OPTIONAL.}
  \item{ylabAngle}{Rotation angle of y-axis labels. DEFAULT = 0. OPTIONAL.}
  \item{ylabhjust}{Horizontal adjustment of y-axis labels. DEFAULT = 0.5.
  OPTIONAL.}
  \item{ylabvjust}{Vertical adjustment of y-axis labels. DEFAULT = 0.5.
  OPTIONAL.}
  \item{axisLabSize}{Size of x- and y-axis labels. DEFAULT = 16. OPTIONAL.}
  \item{title}{Plot title. DEFAULT = ''. OPTIONAL.}
  \item{subtitle}{Plot subtitle. DEFAULT = ''. OPTIONAL.}
  \item{caption}{Plot caption. DEFAULT = ''. OPTIONAL.}
  \item{titleLabSize}{Size of plot title. DEFAULT = 16. OPTIONAL.}
  \item{subtitleLabSize}{Size of plot subtitle. DEFAULT = 12. OPTIONAL.}
  \item{captionLabSize}{Size of plot caption. DEFAULT = 12. OPTIONAL.}
  \item{hline}{Draw one or more horizontal lines passing through this/these
  values on y-axis. For single values, only a single numerical value is
  necessary. For multiple lines, pass these as a vector, e.g., c(60,90).
  DEFAULT = c(0). OPTIONAL.}
  \item{hlineType}{Line type for hline ('blank', 'solid', 'dashed', 'dotted',
  'dotdash', 'longdash', 'twodash'). DEFAULT = 'longdash'. OPTIONAL.}
  \item{hlineCol}{Colour of hline. DEFAULT = 'black'. OPTIONAL.}
  \item{hlineWidth}{Width of hline. DEFAULT = 0.4. OPTIONAL.}
  \item{vline}{Draw one or more vertical lines passing through this/these
  values on x-axis. For single values, only a single numerical value is
  necessary. For multiple lines, pass these as a vector, e.g., c(60,90).
  DEFAULT = NULL. OPTIONAL.}
  \item{vlineType}{Line type for vline ('blank', 'solid', 'dashed', 'dotted',
  'dotdash', 'longdash', 'twodash'). DEFAULT = 'longdash'. OPTIONAL.}
  \item{vlineCol}{Colour of vline. DEFAULT = 'black'. OPTIONAL.}
  \item{vlineWidth}{Width of vline. DEFAULT = 0.4. OPTIONAL.}
  \item{gridlines.major}{Draw major gridlines? (TRUE/FALSE). DEFAULT = TRUE.
  OPTIONAL.}
  \item{gridlines.minor}{Draw minor gridlines? (TRUE/FALSE). DEFAULT = TRUE.
  OPTIONAL.}
  \item{borderWidth}{Width of the border on the x and y axes. DEFAULT = 0.8.
  OPTIONAL.}
  \item{borderColour}{Colour of the border on the x and y axes. DEFAULT =
  'black'. OPTIONAL.}
  \item{returnPlot}{Return the plot object? (TRUE/FALSE). DEFAULT = TRUE.
  OPTIONAL.}
}

\details{Plot the component loadings for selected principal components / eigenvectors and label variables driving variation along these.}

\value{
A \code{\link{ggplot2}} object.
}

\author{
Kevin Blighe <kevin@clinicalbioinformatics.co.uk>, Myles Lewis <myles.lewis@qmul.ac.uk>
}

\examples{
  options(scipen=10)
  options(digits=6)

  col <- 20
  row <- 20000
  mat1 <- matrix(
    rexp(col*row, rate = 0.1),
    ncol = col)
  rownames(mat1) <- paste0('gene', 1:nrow(mat1))
  colnames(mat1) <- paste0('sample', 1:ncol(mat1))

  mat2 <- matrix(
    rexp(col*row, rate = 0.1),
    ncol = col)
  rownames(mat2) <- paste0('gene', 1:nrow(mat2))
  colnames(mat2) <- paste0('sample', (ncol(mat1)+1):(ncol(mat1)+ncol(mat2)))

  mat <- cbind(mat1, mat2)

  metadata <- data.frame(row.names = colnames(mat))
  metadata$Group <- rep(NA, ncol(mat))
  metadata$Group[seq(1,40,2)] <- 'A'
  metadata$Group[seq(2,40,2)] <- 'B'
  metadata$CRP <- sample.int(100, size=ncol(mat), replace=TRUE)
  metadata$ESR <- sample.int(100, size=ncol(mat), replace=TRUE)

  p <- pca(mat, metadata = metadata, removeVar = 0.1)

  plotloadings(p, drawConnectors=TRUE)
}
