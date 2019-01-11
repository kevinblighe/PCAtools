\name{biplot}

\alias{biplot}

\title{biplot}

\description{Draw a bi-plot, comparing 2 selected principal components / eigenvectors.}

\usage{
  biplot(pcaobj,
  x = 'PC1',
  y = 'PC2',
  colby = NULL,
  colkey = NULL,
  singlecol = NULL,
  shape = NULL,
  shapekey = NULL,
  pointSize = 3.0,
  legendPosition = 'none',
  legendLabSize = 12,
  legendIconSize = 5.0,
  xlim = NULL,
  ylim = NULL,
  lab = TRUE,
  labSize = 3.0,
  labhjust = 1.5,
  labvjust = 0,
  selectLab = NULL,
  drawConnectors = TRUE,
  widthConnectors = 0.5,
  colConnectors = 'grey50',
  xlab = paste0(x, ', ',
    round(pcaobj$variance[x], digits=2),
    '\% variation'),
  xlabAngle = 0,
  xlabhjust = 0.5,
  xlabvjust = 0.5,
  ylab = paste0(y, ', ',
    round(pcaobj$variance[y], digits=2),
    '\% variation'),
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
  hline = NULL,
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
  \item{x}{A principal component to plot on x-axis. All principal component
  names are stored in pcaobj$label. DEFAULT = 'PC1'. REQUIRED.}
  \item{y}{A principal component to plot on y-axis. All principal component
  names are stored in pcaobj$label. DEFAULT = 'PC2'. REQUIRED.}
  \item{colby}{If NULL, all points will be coloured differently. If not NULL,
  value is assumed to be a column name in pcaobj$metadata relating to some
  grouping/categorical variable. DEFAULT = NULL. OPTIONAL.}
  \item{colkey}{Vector of name-value pairs relating to value passed to 'col',
  e.g., c(A='forestgreen', B='gold'). DEFAULT = NULL. OPTIONAL.}
  \item{singlecol}{If specified, all points will be shaded by this colour.
  Overrides 'col'. DEFAULT = NULL. OPTIONAL.}
  \item{shape}{If NULL, all points will be have the same shape. If not NULL,
  value is assumed to be a column name in pcaobj$metadata relating to some
  grouping/categorical variable. DEFAULT = NULL. OPTIONAL.}
  \item{shapekey}{Vector of name-value pairs relating to value passed to
  'shape', e.g., c(A=10, B=21). DEFAULT = NULL. OPTIONAL.}
  \item{pointSize}{Size of plotted points. DEFAULT = 3.0. OPTIONAL.}
  \item{legendPosition}{Position of legend ('top', 'bottom', 'left', 'right',
  'none'). DEFAULT = 'none'. OPTIONAL.}
  \item{legendLabSize}{Size of plot legend text. DEFAULT = 10. OPTIONAL.}
  \item{legendIconSize}{Size of plot legend icons / symbols. DEFAULT = 3.0.
  OPTIONAL.}
  \item{xlim}{Limits of the x-axis. DEFAULT = NULL. OPTIONAL.}
  \item{ylim}{Limits of the y-axis. DEFAULT = NULL. OPTIONAL.}
  \item{lab}{Logical, indicating whether or not to label the points in the
  plot space. Labels will be taken as the original colnames of the input
  object, usually sample IDs. DEFAULT = TRUE. OPTIONAL.}
  \item{labSize}{Size of labels. DEFAULT = 3.0. OPTIONAL.}
  \item{labhjust}{Horizontal adjustment of label. DEFAULT = 1.5. OPTIONAL.}
  \item{labvjust}{Vertical adjustment of label. DEFAULT = 0. OPTIONAL.}
  \item{selectLab}{A vector containing a subset of lab to plot. DEFAULT =
  NULL. OPTIONAL.}
  \item{drawConnectors}{Logical, indicating whether or not to connect plot
  labels to their corresponding points by line connectors. DEFAULT = TRUE.
  OPTIONAL.}
  \item{widthConnectors}{Line width of connectors. DEFAULT = 0.5. OPTIONAL.}
  \item{colConnectors}{Line colour of connectors. DEFAULT = 'grey50'. OPTIONAL.}
  \item{xlab}{Label for x-axis. DEFAULT =
  paste0(x, ', ', round(pcaobj$variance[x], digits=2), '\% variation').
  OPTIONAL.}
  \item{xlabAngle}{Rotation angle of x-axis labels. DEFAULT = 0. OPTIONAL.}
  \item{xlabhjust}{Horizontal adjustment of x-axis labels. DEFAULT = 0.5. OPTIONAL.}
  \item{xlabvjust}{Vertical adjustment of x-axis labels. DEFAULT = 0.5.
  OPTIONAL.}
  \item{ylab}{Label for y-axis. DEFAULT = paste0(y, ', ',
  round(pcaobj$variance[y], digits=2), '\% variation'). OPTIONAL.}
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
  DEFAULT = NULL. OPTIONAL.}
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
  \item{gridlines.major}{Logical, indicating whether or not to draw major
  gridlines. DEFAULT = TRUE. OPTIONAL.}
  \item{gridlines.minor}{Logical, indicating whether or not to draw minor
  gridlines. DEFAULT = TRUE. OPTIONAL.}
  \item{borderWidth}{Width of the border on the x and y axes. DEFAULT = 0.8.
  OPTIONAL.}
  \item{borderColour}{Colour of the border on the x and y axes. DEFAULT =
  'black'. OPTIONAL.}
  \item{returnPlot}{Logical, indicating whether or not to return the plot
  object. DEFAULT = TRUE. OPTIONAL.}
}

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

  biplot(p)

  biplot(p, colby = 'Group', shape = 'Group')

  biplot(p, colby = 'Group', colkey = c(A = 'forestgreen', B = 'gold'),
    legendPosition = 'right')

  biplot(p, colby = 'Group', colkey = c(A='forestgreen', B='gold'),
    shape = 'Group', shapekey = c(A=10, B=21), legendPosition = 'bottom')

}
