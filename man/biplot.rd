\name{biplot}

\alias{biplot}

\title{biplot}

\description{Draw a bi-plot, comparing 2 selected principal components / eigenvectors.}

\usage{biplot(
  pcaobj,
  x = 'PC1',
  y = 'PC2',
  showLoadings = FALSE,
  ntopLoadings = 5,
  showLoadingsNames = if (showLoadings) TRUE else FALSE,
  colLoadingsNames = 'black',
  sizeLoadingsNames = 3,
  boxedLoadingsNames = TRUE,
  fillBoxedLoadings = alpha('white', 1/4),
  drawConnectorsLoadings = TRUE,
  widthConnectorsLoadings = 0.5,
  colConnectorsLoadings = 'grey50',
  lengthLoadingsArrowsFactor = 1.5,
  colLoadingsArrows = 'black',
  widthLoadingsArrows = 0.5,
  alphaLoadingsArrow = 1.0,
  colby = NULL,
  colkey = NULL,
  colLegendTitle = if (!is.null(colby)) colby else NULL,
  singlecol = NULL,
  shape = NULL,
  shapekey = NULL,
  shapeLegendTitle = if (!is.null(shape)) shape else NULL,
  pointSize = 3.0,
  legendPosition = 'none',
  legendLabSize = 12,
  legendTitleSize = 14,
  legendIconSize = 5.0,
  encircle = FALSE,
  encircleFill = TRUE,
  encircleFillKey = NULL,
  encircleAlpha = 1/4,
  encircleLineSize = 0.25,
  encircleLineCol = NULL,
  ellipse = FALSE,
  ellipseConf = 0.95,
  ellipseFill = TRUE,
  ellipseFillKey = NULL,
  ellipseAlpha = 1/4,
  ellipseLineSize = 0.25,
  ellipseLineCol = NULL,
  xlim = if(showLoadings) c(min(pcaobj$rotated[,x]) - 5, max(pcaobj$rotated[,x]) + 5) else c(min(pcaobj$rotated[,x]) - 1, max(pcaobj$rotated[,x]) + 1),
  ylim = if(showLoadings) c(min(pcaobj$rotated[,y]) - 5, max(pcaobj$rotated[,y]) + 5) else c(min(pcaobj$rotated[,y]) - 1, max(pcaobj$rotated[,y]) + 1),
  lab = rownames(pcaobj$metadata),
  labSize = 3.0,
  labhjust = 1.5,
  labvjust = 0,
  boxedLabels = FALSE,
  selectLab = NULL,
  drawConnectors = TRUE,
  widthConnectors = 0.5,
  colConnectors = 'grey50',
  xlab = paste0(x, ', ', round(pcaobj$variance[x], digits = 2), '% variation'),
  xlabAngle = 0,
  xlabhjust = 0.5,
  xlabvjust = 0.5,
  ylab = paste0(y, ', ', round(pcaobj$variance[y], digits = 2), '% variation'),
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
  \item{showLoadings}{Logical, indicating whether or not to overlay
    variable loadings. DEFAULT = FALSE. OPTIONAL.}
  \item{ntopLoadings}{If showLoadings == TRUE, select this many variables
    based on absolute ordered variable loading for each PC in the biplot.
    As a result of looking across 2 PCs, it can occur whereby greater than
    this number are actually displayed. DEFAULT = 5. OPTIONAL.}
  \item{showLoadingsNames}{Logical, indicating to show variable loadings names
    or not. DEFAULT = if (showLoadings) TRUE else FALSE. OPTIONAL.}
  \item{colLoadingsNames}{If 'showLoadings == TRUE', colour of text labels.
    DEFAULT = 'black'. OPTIONAL.}
  \item{sizeLoadingsNames}{If 'showLoadings == TRUE', size of text labels.
    DEFAULT = 3. OPTIONAL.}
  \item{boxedLoadingsNames}{Logical, if 'showLoadings == TRUE', draw text
    labels in boxes. DEFAULT = TRUE. OPTIONAL.}
  \item{fillBoxedLoadings}{When 'boxedLoadingsNames == TRUE', this controls
    the background fill of the boxes. To control both the fill and
    transparency, user can specify a value of the form
    'alpha(<colour>, <alpha>)'. DEFAULT = alpha('white', 1/4). OPTIONAL.}
  \item{drawConnectorsLoadings}{If 'showLoadings == TRUE', draw line connectors
    to the variable loadings arrows in order to fit more labels in the plot
    space. DEFAULT = TRUE. OPTIONAL.}
  \item{widthConnectorsLoadings}{If 'showLoadings == TRUE', width of the line
    connectors drawn to the variable loadings arrows. DEFAULT = 0.5. OPTIONAL.}
  \item{colConnectorsLoadings}{If 'showLoadings == TRUE', colour of the line
    connectors drawn to the variable loadings arrows. DEFAULT = 'grey50'.
    OPTIONAL.}
  \item{lengthLoadingsArrowsFactor}{If 'showLoadings == TRUE', multiply the
    internally-determined length of the variable loadings arrows by this
    factor. DEFAULT = 1.5. OPTIONAL.}
  \item{colLoadingsArrows}{If showLoadings == TRUE, colour of the variable
    loadings arrows. DEFAULT = 'black'. OPTIONAL.}
  \item{widthLoadingsArrows}{If showLoadings == TRUE, width of the variable
    loadings arrows. DEFAULT = 0.5. OPTIONAL.}
  \item{alphaLoadingsArrow}{If showLoadings == TRUE, colour transparency of
    the variable loadings arrows. DEFAULT = 1.0. OPTIONAL.}
  \item{colby}{If NULL, all points will be coloured differently. If not NULL,
    value is assumed to be a column name in pcaobj$metadata relating to some
    grouping/categorical variable. DEFAULT = NULL. OPTIONAL.}
  \item{colkey}{Vector of name-value pairs relating to value passed to 'col',
    e.g., c(A='forestgreen', B='gold'). DEFAULT = NULL. OPTIONAL.}
  \item{colLegendTitle}{Title of the legend for the variable specified
    by 'colby'. DEFAULT = if (!is.null(colby)) colby else NULL. OPTIONAL.}
  \item{singlecol}{If specified, all points will be shaded by this colour.
    Overrides 'col'. DEFAULT = NULL. OPTIONAL.}
  \item{shape}{If NULL, all points will be have the same shape. If not NULL,
    value is assumed to be a column name in pcaobj$metadata relating to some
    grouping/categorical variable. DEFAULT = NULL. OPTIONAL.}
  \item{shapekey}{Vector of name-value pairs relating to value passed to
    'shape', e.g., c(A=10, B=21). DEFAULT = NULL. OPTIONAL.}
  \item{shapeLegendTitle}{Title of the legend for the variable specified
    by 'shape'. DEFAULT = if (!is.null(shape)) shape else NULL. OPTIONAL.}
  \item{pointSize}{Size of plotted points. DEFAULT = 3.0. OPTIONAL.}
  \item{legendPosition}{Position of legend ('top', 'bottom', 'left', 'right',
    'none'). DEFAULT = 'none'. OPTIONAL.}
  \item{legendLabSize}{Size of plot legend text. DEFAULT = 12. OPTIONAL.}
  \item{legendTitleSize}{Size of plot legend title text. DEFAULT = 14. OPTIONAL.}
  \item{legendIconSize}{Size of plot legend icons / symbols. DEFAULT = 5.0.
    OPTIONAL.}
  \item{encircle}{Logical, indicating whether to draw a polygon around
    the groups specified by 'colby'. DEFAULT = FALSE. OPTIONAL.}
  \item{encircleFill}{Logical, if 'encircle == TRUE', this determines
    whether to fill the encircled region or not. DEFAULT = TRUE. OPTIONAL.}
  \item{encircleFillKey}{Vector of name-value pairs relating to value passed to
    'encircleFill', e.g., c(A='forestgreen', B='gold'). If NULL, the fill
    is controlled by whatever has already been used for 'colby' / 'colkey'.
    DEFAULT = NULL. OPTIONAL.}
  \item{encircleAlpha}{Alpha for purposes of controlling colour transparency of
    the encircled region. Used when 'encircle == TRUE'. DEFAULT = 1/4.
    OPTIONAL.}
  \item{encircleLineSize}{Line width of the encircled line when
    'encircle == TRUE'. DEFAULT = 0.25. OPTIONAL.}
  \item{encircleLineCol}{Colour of the encircled line when
    'encircle == TRUE'. DEFAULT = NULL. OPTIONAL.}
  \item{ellipse}{Logical, indicating whether to draw a stat ellipse around
    the groups specified by 'colby'. DEFAULT = FALSE. OPTIONAL.}
  \item{ellipseConf}{Confidence intervals of the stat ellipses when
    ellipse == TRUE. DEFAULT = 0.95. OPTIONAL.}
  \item{ellipseFill}{Logical, if 'ellipse == TRUE', this determines
    whether to fill the region or not. DEFAULT = TRUE. OPTIONAL.}
  \item{ellipseFillKey}{Vector of name-value pairs relating to value passed to
    'ellipseFill', e.g., c(A='forestgreen', B='gold'). If NULL, the fill
    is controlled by whatever has already been used for 'colby' / 'colkey'.
    DEFAULT = NULL. OPTIONAL.}
  \item{ellipseAlpha}{Alpha for purposes of controlling colour transparency of
    the ellipse region. Used when 'ellipse == TRUE'. DEFAULT = 1/4. OPTIONAL.}
  \item{ellipseLineSize}{Line width of the ellipse line when 'ellipse == TRUE'.
    DEFAULT = 0.25. OPTIONAL.}
  \item{ellipseLineCol}{Colour of the ellipse line when 'ellipse == TRUE'.
    DEFAULT = NULL. OPTIONAL.}
  \item{xlim}{Limits of the x-axis. DEFAULT = if(showLoadings)
    c(min(pcaobj$rotated[,x]) - 5, max(pcaobj$rotated[,x]) + 5)
    else c(min(pcaobj$rotated[,x]) - 1, max(pcaobj$rotated[,x]) + 1).
    OPTIONAL.}
  \item{ylim}{Limits of the y-axis. DEFAULT = if(showLoadings)
    c(min(pcaobj$rotated[,y]) - 5, max(pcaobj$rotated[,y]) + 5)
    else c(min(pcaobj$rotated[,y]) - 1, max(pcaobj$rotated[,y]) + 1).
    OPTIONAL.}
  \item{lab}{A vector containing labels to add to the plot. 
    DEFAULT = rownames(pcaobj$metadata). OPTIONAL.}
  \item{labSize}{Size of labels. DEFAULT = 3.0. OPTIONAL.}
  \item{labhjust}{Horizontal adjustment of label. DEFAULT = 1.5. OPTIONAL.}
  \item{labvjust}{Vertical adjustment of label. DEFAULT = 0. OPTIONAL.}
  \item{boxedLabels}{Logical, draw text labels in boxes. DEFAULT = FALSE.
    OPTIONAL.}
  \item{selectLab}{A vector containing a subset of lab to plot. DEFAULT =
    NULL. OPTIONAL.}
  \item{drawConnectors}{Logical, indicating whether or not to connect plot
    labels to their corresponding points by line connectors. DEFAULT = TRUE.
    OPTIONAL.}
  \item{widthConnectors}{Line width of connectors. DEFAULT = 0.5. OPTIONAL.}
  \item{colConnectors}{Line colour of connectors. DEFAULT = 'grey50'.
    OPTIONAL.}
  \item{xlab}{Label for x-axis. DEFAULT = paste0(x, ', ',
    round(pcaobj$variance[x], digits = 2), '% variation'). OPTIONAL.}
  \item{xlabAngle}{Rotation angle of x-axis labels. DEFAULT = 0. OPTIONAL.}
  \item{xlabhjust}{Horizontal adjustment of x-axis labels. DEFAULT = 0.5.
    OPTIONAL.}
  \item{xlabvjust}{Vertical adjustment of x-axis labels. DEFAULT = 0.5.
    OPTIONAL.}
  \item{ylab}{Label for y-axis. DEFAULT = paste0(y, ', ',
    round(pcaobj$variance[y], digits = 2), '% variation'). OPTIONAL.}
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

\value{A \code{\link{ggplot2}} object.}

\author{Kevin Blighe <kevin@clinicalbioinformatics.co.uk>}

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
