\name{screeplot}

\alias{screeplot}

\title{screeplot}

\description{Draw a SCREE plot, showing the distribution of explained variance across all or select principal components / eigenvectors.}

\usage{
  screeplot(pcaobj,
  components = getComponents(pcaobj),
  xlim = NULL,
  ylim = c(0, 100),
  xlab = 'Principal component',
  xlabAngle = 90,
  xlabhjust = 0.5,
  xlabvjust = 0.5,
  ylab = 'Explained variation (\%)',
  ylabAngle = 0,
  ylabhjust = 0.5,
  ylabvjust = 0.5,
  axisLabSize = 16,
  title = 'SCREE plot',
  subtitle = '',
  caption = '',
  titleLabSize = 16,
  subtitleLabSize = 12,
  captionLabSize = 12,
  colBar = 'dodgerblue',
  drawCumulativeSumLine = TRUE,
  colCumulativeSumLine = 'red2',
  sizeCumulativeSumLine = 1.5,
  drawCumulativeSumPoints = TRUE,
  colCumulativeSumPoints = 'red2',
  sizeCumulativeSumPoints = 2.0,
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
  \item{components}{The principal components to be included in the plot.
  DEFAULT = getComponents(pcaobj). OPTIONAL.}
  \item{xlim}{Limits of the x-axis. DEFAULT = NULL. OPTIONAL.}
  \item{ylim}{Limits of the y-axis. DEFAULT = c(0, 100). OPTIONAL.}
  \item{xlab}{Label for x-axis. DEFAULT = 'Principal component'. OPTIONAL.}
  \item{xlabAngle}{Rotation angle of x-axis labels. DEFAULT = 90. OPTIONAL.}
  \item{xlabhjust}{Horizontal adjustment of x-axis labels. DEFAULT = 0.5.
  OPTIONAL.}
  \item{xlabvjust}{Vertical adjustment of x-axis labels. DEFAULT = 0.5.
  OPTIONAL.}
  \item{ylab}{Label for y-axis. DEFAULT = 'Explained variation (\%)'. OPTIONAL.}
  \item{ylabAngle}{Rotation angle of y-axis labels. DEFAULT = 0. OPTIONAL.}
  \item{ylabhjust}{Horizontal adjustment of y-axis labels. DEFAULT = 0.5.
  OPTIONAL.}
  \item{ylabvjust}{Vertical adjustment of y-axis labels. DEFAULT = 0.5.
  OPTIONAL.}
  \item{axisLabSize}{Size of x- and y-axis labels. DEFAULT = 16. OPTIONAL.}
  \item{title}{Plot title. DEFAULT = 'SCREE plot'. OPTIONAL.}
  \item{subtitle}{Plot subtitle. DEFAULT = ''. OPTIONAL.}
  \item{caption}{Plot caption. DEFAULT = ''. OPTIONAL.}
  \item{titleLabSize}{Size of plot title. DEFAULT = 16. OPTIONAL.}
  \item{subtitleLabSize}{Size of plot subtitle. DEFAULT = 12. OPTIONAL.}
  \item{captionLabSize}{Size of plot caption. DEFAULT = 12. OPTIONAL.}
  \item{colBar}{DEFAULT = 'dodgerblue'. OPTIONAL.}
  \item{drawCumulativeSumLine}{Logical, indicating whether or not to overlay
  plot with a cumulative explained variance line. DEFAULT = TRUE. OPTIONAL.}
  \item{colCumulativeSumLine}{Colour of cumulative explained variance line.
  DEFAULT = 'red2'. OPTIONAL.}
  \item{sizeCumulativeSumLine}{Size of cumulative explained variance line.
  DEFAULT = 1.5. OPTIONAL.}
  \item{drawCumulativeSumPoints}{Logical, indicating whether or not to draw
  the cumulative explained variance points. DEFAULT = TRUE. OPTIONAL.}
  \item{colCumulativeSumPoints}{Colour of cumulative explained variance
  points. DEFAULT = 'red2'. OPTIONAL.}
  \item{sizeCumulativeSumPoints}{Size of cumulative explained variance points.
  DEFAULT = 2.0. OPTIONAL.}
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
Kevin Blighe <kevin@clinicalbioinformatics.co.uk>
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

  screeplot(p)

  screeplot(p, hline = 80)
}
