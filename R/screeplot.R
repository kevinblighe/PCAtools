#' Draw a SCREE plot, showing the distribution of explained variance across all or select principal components / eigenvectors.
#'
#' @param pcaobj Object of class 'pca' created by pca().
#' @param components The principal components to be included in the plot.
#' @param xlim Limits of the x-axis.
#' @param ylim Limits of the y-axis.
#' @param xlab Label for x-axis.
#' @param xlabAngle Rotation angle of x-axis labels.
#' @param xlabhjust Horizontal adjustment of x-axis labels.
#' @param xlabvjust Vertical adjustment of x-axis labels.
#' @param ylab Label for y-axis.
#' @param ylabAngle Rotation angle of y-axis labels.
#' @param ylabhjust Horizontal adjustment of y-axis labels.
#' @param ylabvjust Vertical adjustment of y-axis labels.
#' @param axisLabSize Size of x- and y-axis labels.
#' @param title Plot title.
#' @param subtitle Plot subtitle.
#' @param caption Plot caption.
#' @param titleLabSize Size of plot title.
#' @param subtitleLabSize Size of plot subtitle.
#' @param captionLabSize Size of plot caption.
#' @param colBar Colour of the vertical bars.
#' @param drawCumulativeSumLine Logical, indicating whether or not to overlay
#'   plot with a cumulative explained variance line.
#' @param colCumulativeSumLine Colour of cumulative explained variance line.
#' @param sizeCumulativeSumLine Size of cumulative explained variance line.
#' @param drawCumulativeSumPoints Logical, indicating whether or not to draw
#'   the cumulative explained variance points.
#' @param colCumulativeSumPoints Colour of cumulative explained variance
#'   points.
#' @param sizeCumulativeSumPoints Size of cumulative explained variance points.
#' @param hline Draw one or more horizontal lines passing through this/these
#'   values on y-axis. For single values, only a single numerical value is
#'   necessary. For multiple lines, pass these as a vector, e.g., c(60,90).
#' @param hlineType Line type for hline ('blank', 'solid', 'dashed', 'dotted',
#'   'dotdash', 'longdash', 'twodash').
#' @param hlineCol Colour of hline.
#' @param hlineWidth Width of hline.
#' @param vline Draw one or more vertical lines passing through this/these
#'   values on x-axis. For single values, only a single numerical value is
#'   necessary. For multiple lines, pass these as a vector, e.g., c(60,90).
#' @param vlineType Line type for vline ('blank', 'solid', 'dashed', 'dotted',
#'   'dotdash', 'longdash', 'twodash').
#' @param vlineCol Colour of vline.
#' @param vlineWidth Width of vline.
#' @param gridlines.major Logical, indicating whether or not to draw major
#'   gridlines.
#' @param gridlines.minor Logical, indicating whether or not to draw minor
#'   gridlines.
#' @param borderWidth Width of the border on the x and y axes.
#' @param borderColour Colour of the border on the x and y axes.
#' @param returnPlot Logical, indicating whether or not to return the plot
#'   object.
#'
#' @details Draw a SCREE plot, showing the distribution of explained variance across all or select principal components / eigenvectors.
#'
#' @return A \code{\link{ggplot2}} object.
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
#'   screeplot(p)
#'
#'   screeplot(p, hline = 80)
#'
#' @import ggplot2
#'
#' @export
screeplot <- function(
  pcaobj,
  components = getComponents(pcaobj),
  xlim = NULL,
  ylim = c(0, 100),
  xlab = 'Principal component',
  xlabAngle = 90,
  xlabhjust = 0.5,
  xlabvjust = 0.5,
  ylab = 'Explained variation (%)',
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
{
  PC <- Variance <- NULL

  # create a plot object (2-col df) of PC names and
  # explained variance of each
  plotobj <- data.frame(components, pcaobj$variance[components])
  colnames(plotobj) <- c('PC', 'Variance')
  plotobj$PC <- factor(plotobj$PC,
    levels=plotobj$PC[seq_along(plotobj$PC)])

  # create a base theme that will later be modified
  th <- theme_bw(base_size=24) +

    theme(
      legend.background=element_rect(),
      plot.title=element_text(angle = 0, size = titleLabSize,
        face = 'bold', vjust = 1),
      plot.subtitle=element_text(angle = 0, size = subtitleLabSize,
        face = 'plain', vjust = 0.5),
      plot.caption=element_text(angle = 0, size = captionLabSize,
        face = 'plain', vjust = 0.5),

      axis.text.x=element_text(angle = xlabAngle, size = axisLabSize,
        hjust = xlabhjust, vjust = xlabvjust),
      axis.text.y=element_text(angle = ylabAngle, size = axisLabSize,
        hjust = ylabhjust, vjust = ylabvjust),
      axis.title=element_text(size = axisLabSize),

      legend.position = 'none')

  # create the plot as geom_bar
  plot <- ggplot(plotobj, aes(x = PC)) + th +
    geom_bar(aes(y = Variance), fill = colBar,
      width = 0.8, stat = 'identity', show.legend = FALSE)

  # add cumulative sum line as geom_line
  if (drawCumulativeSumLine == TRUE) {
    plot <- plot + geom_line(aes(y = cumsum(Variance), group = 1),
      colour = colCumulativeSumLine,
      size = sizeCumulativeSumLine,
      show.legend = FALSE)
  }

  # add cumulative sum points as geom_point
  if (drawCumulativeSumPoints == TRUE) {
    plot <- plot + geom_point(aes(y = cumsum(Variance)),
      colour = colCumulativeSumPoints,
      size = sizeCumulativeSumPoints,
      show.legend = FALSE)
  }

  # add elements to the plot for xy labeling and axis limits
  plot <- plot + xlab(xlab) + ylab(ylab)
  if (!is.null(xlim)) {
    plot <- plot + xlim(xlim[1], xlim[2])
  }
  if (!is.null(ylim)) {
    plot <- plot + ylim(ylim[1], ylim[2])
  }

  # add elements to the plot for title, subtitle, caption
  plot <- plot + labs(title = title, 
    subtitle = subtitle, caption = caption)

  # add elements to the plot for vlines and hlines
  if (!is.null(vline)) {
    plot <- plot + geom_vline(xintercept = vline,
      linetype = vlineType,
      colour = vlineCol,
      size = vlineWidth)
  }
  if (!is.null(hline)) {
    plot <- plot + geom_hline(yintercept = hline,
      linetype = hlineType,
      colour = hlineCol,
      size = hlineWidth)
  }

  # border around plot
  plot <- plot +
    theme(panel.border = element_rect(
      colour = borderColour,
      fill = NA,
      size = borderWidth))

  # gridlines
  if (gridlines.major == TRUE) {
    plot <- plot + theme(panel.grid.major = element_line())
  } else {
    plot <- plot + theme(panel.grid.major = element_blank())
  }
  if (gridlines.minor == TRUE) {
    plot <- plot + theme(panel.grid.minor = element_line())
  } else {
    plot <- plot + theme(panel.grid.minor = element_blank())
  }

  # return plot?
  if (returnPlot == TRUE) {
    return(plot)
  } else if (returnPlot == FALSE) {
    plot
  }
}
