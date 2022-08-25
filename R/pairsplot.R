#' Draw multiple bi-plots.
#'
#' @param pcaobj Object of class 'pca' created by pca().
#' @param components The principal components to be included in the plot. These
#'   will be compared in a pairwise fashion via multiple calls to biplot().
#' @param triangle Logical, indicating whether or not to draw the plots in the
#'   upper panel in a triangular arrangement? Principal component names will be
#'   labeled along the diagonal.
#' @param trianglelabSize Size of p rincipal component label (when triangle =
#'   TRUE).
#' @param plotaxes Logical, indicating whether or not to draw the axis tick,
#'   labels, and titles.
#' @param margingaps The margins between plots in the plot space. Takes the form
#'   of a 'unit()' variable.
#' @param ncol If triangle = FALSE, the number of columns in the final merged
#'   plot.
#' @param nrow If triangle = FALSE, the number of rows in the final merged
#'   plot.
#' @param x A principal component to plot on x-axis. All principal component
#'   names are stored in pcaobj$label.
#' @param y A principal component to plot on y-axis. All principal component
#'   names are stored in pcaobj$label.
#' @param colby If NULL, all points will be coloured differently. If not NULL,
#'   value is assumed to be a column name in pcaobj$metadata relating to some
#'   grouping/categorical variable.
#' @param colkey Vector of name-value pairs relating to value passed to 'col',
#'   e.g., c(A='forestgreen', B='gold').
#' @param singlecol If specified, all points will be shaded by this colour.
#'   Overrides 'col'.
#' @param shape If NULL, all points will be have the same shape. If not NULL,
#'   value is assumed to be a column name in pcaobj$metadata relating to some
#'   grouping/categorical variable.
#' @param shapekey Vector of name-value pairs relating to value passed to
#'   'shape', e.g., c(A=10, B=21).
#' @param pointSize Size of plotted points.
#' @param shared.legend Logical, indicating whether to draw a shared legend for the bi-plots.
#' @param legendPosition Position of legend ('top', 'bottom', 'left', 'right',
#'   'none').
#' @param legendLabSize Size of plot legend text.
#' @param legendIconSize Size of plot legend icons / symbols.
#' @param xlim Limits of the x-axis.
#' @param ylim Limits of the y-axis.
#' @param lab A vector containing labels to add to the plot.
#' @param labSize Size of labels.
#' @param selectLab A vector containing a subset of lab to plot.
#' @param drawConnectors Logical, indicating whether or not to connect plot
#'   labels to their corresponding points by line connectors.
#' @param widthConnectors Line width of connectors.
#' @param colConnectors Line colour of connectors.
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
#' @param titleLabSize Size of plot title.
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
#' @details Draw multiple bi-plots.
#'
#' @return A \code{\link{cowplot}} object.
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
#'   pairsplot(p, triangle = TRUE)
#'
#' @import ggplot2
#' @import cowplot
#'
#' @export
pairsplot <- function(
  pcaobj,
  components = getComponents(pcaobj, seq_len(5)),
  triangle = TRUE,
  trianglelabSize = 18,
  plotaxes = TRUE,
  margingaps = unit(c(0.1, 0.1, 0.1, 0.1), 'cm'),
  ncol = NULL,
  nrow = NULL,
  # other biplot() params:
  x = NULL,
  y = NULL,
  colby = NULL,
  colkey = NULL,
  singlecol = NULL,
  shape = NULL,
  shapekey = NULL,
  pointSize = 1.0,
  shared.legend = FALSE,
  legendPosition = 'none',
  legendLabSize = 6,
  legendTitleSize = 14,
  legendIconSize = 1.5,
  xlim = NULL,
  ylim = NULL,
  lab = NULL,
  labSize = 1.5,
  selectLab = NULL,
  drawConnectors = FALSE,
  widthConnectors = 0.5,
  colConnectors = 'grey50',
  xlab = NULL,
  xlabAngle = 0,
  xlabhjust = 0.5,
  xlabvjust = 0.5,
  ylab = NULL,
  ylabAngle = 0,
  ylabhjust = 0.5,
  ylabvjust = 0.5,
  axisLabSize = 10,
  title = NULL,
  titleLabSize = 32,
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
  # biplots is a list that will be populated with biplot()
  # function return objects and any geom labels
  biplots <- list()

  # counter necessary for layout of objects in plot space
  nplots <- 0
  
  # the shared legend will be on the left
  if (shared.legend) {
    legendPosition <- 'left'
  }
  
  # beginning of the master loop (contains nested loop)
  # biplots will be created on a pairwise basis
  for (i in seq_along(components)) {

    # if the triangular layout is specified, large titles are added
    # as ggdraw() objects
    if (triangle == TRUE) {
      nplots <- nplots + 1 # increment nplots

      biplots[[nplots]] <- ggdraw() +
        draw_label(x = 0.6, y = 0.6,
          paste0(components[i],
            ',\n',
            round(pcaobj$variance[components[i]], digits=2), '%'),
          fontface = 'bold',
          size = trianglelabSize)

      names(biplots)[nplots] <- 'Label'
    }

    # nested loop
    for (j in seq_along(components)) {

      # if statement that prevents 'self' and duplicate biplots
      if (i != j && i < j) {
        nplots <- nplots + 1 # increment nplots

        x <- components[i]
        y <- components[j]

        # call the biplot() function
        biplots[[nplots]] <- biplot(pcaobj,
          x = x,
          y = y,
          # other biplot() params:
          colby = colby,
          colkey = colkey,
          singlecol = singlecol,
          shape = shape,
          shapekey = shapekey,
          selectLab = selectLab,
          xlim = xlim,
          ylim = ylim,
          lab = lab,
          xlab = paste0(x, ', ',
            round(pcaobj$variance[x], digits=2),
            '%'),
          xlabAngle = xlabAngle,
          xlabhjust = xlabhjust,
          xlabvjust = xlabvjust,
          ylab = paste0(y, ', ',
            round(pcaobj$variance[y], digits=2),
            '%'),
          ylabAngle = ylabAngle,
          ylabhjust = ylabhjust,
          ylabvjust = ylabvjust,
          axisLabSize = axisLabSize,
          pointSize = pointSize,
          labSize = labSize,
          legendPosition = legendPosition,
          legendLabSize = legendLabSize,
          legendTitleSize = legendTitleSize,
          legendIconSize = legendIconSize,
          drawConnectors = drawConnectors,
          widthConnectors = widthConnectors,
          colConnectors = colConnectors,
          hline = hline,
          hlineType = hlineType,
          hlineCol = hlineCol,
          hlineWidth = hlineWidth,
          vline = vline,
          vlineType = vlineType,
          vlineCol = vlineCol,
          vlineWidth = vlineWidth,
          gridlines.major = gridlines.major,
          gridlines.minor = gridlines.minor,
          borderWidth = borderWidth,
          borderColour = borderColour,
          returnPlot = returnPlot)

       # assign list name to plot, e.g. 'PC1 Vs PC3'
       names(biplots)[nplots] <- paste(components[i], 'Vs', components[j])
      }
    }
  }

  # specify margin (gaps and titles) that will be added to each plot
  # remove titles, subtitles, caption
  margin <- theme(
    plot.margin = margingaps,
    plot.title = element_blank(),
    plot.subtitle = element_blank(),
    plot.caption = element_blank())

  # save the title as a ggdraw object
  title <- ggdraw() + draw_label(title,
    fontface = 'bold',
    size = titleLabSize)

  # triangular layout?
  # with triangular layout, empty space is filled as empty
  # plots with label = ''
  if (triangle) {
    ncol <- nrow <- length(components)

    # create new list that will store plot objects plus
    # necessary empty plots to form triangular layout
    biplots.final <- list()

    nplots.final <- 0
    nplots.original <- 0

    l <- 0

    for (k in seq_along(components)) {
      # determine number of plots where each PC is on x-axis
      # looks at biplot list names, i.e., 'PC1 Vs PC3', 'PC1 Vs PC4',
      # 'PC1 Vs PC3', 'PC2 Vs PC3', et cetera. Here, numplot=3 for PC1
      numplot <- length(grep(paste0('^', components[k]), names(biplots))) + 1

      # based on numplot, add necessary number of blank plots
      l <- numplot
      while (l < length(components)) {
        nplots.final <- nplots.final + 1
        biplots.final[[nplots.final]] <- ggdraw() +
          draw_label(x = 0.5, y = 0.5, '', fontface = 'bold', size = 32)
        l <- l + 1
      }

      # after adding blank plots, add all original plots
      for (m in seq_len(numplot)) {
        nplots.final <- nplots.final + 1
        nplots.original <- nplots.original + 1
        biplots.final[[nplots.final]] <- biplots[[nplots.original]]
      }
    }

    # remove axis labels and ticks?
    if (!plotaxes) {
      margin <- margin + theme(
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title = element_blank())
    } else if (plotaxes) {
      margin <- margin
    }

    # apply titles / axis changes
    biplots.final <- lapply(biplots.final, '+', margin)
    biplots.final <- lapply(biplots.final, '+', coord_flip())

    # return plot?
    if (returnPlot) {
      if(shared.legend) {
        legend <- get_legend(biplots.final[[2]] + theme(legend.box.margin = margin(0, 0, 0, 12)))
        bigrid <- plot_grid(title, 
                            do.call(plot_grid, c(lapply(biplots.final, "+", theme(legend.position = "none")), ncol = ncol, nrow = nrow)), 
                            ncol = 1, rel_heights = c(0.1, 1))
        return(plot_grid(legend, bigrid, rel_widths = c(0.5, ncol)))
      } else {
        return(plot_grid(title, 
                         do.call(plot_grid, c(biplots.final, ncol = ncol, nrow = nrow)), 
                         ncol = 1, rel_heights = c(0.1, 1)))
      }
    } else if (!returnPlot) {
      if(shared.legend) {
        legend <- get_legend(biplots.final[[2]] + theme(legend.box.margin = margin(0, 0, 0, 12)))
        bigrid <- plot_grid(title, 
                            do.call(plot_grid, c(lapply(biplots.final, "+", theme(legend.position = "none")), ncol = ncol, nrow = nrow)), 
                            ncol = 1, rel_heights = c(0.1, 1))
        plot_grid(legend, bigrid, rel_widths = c(0.5, ncol))
      } else {
        plot_grid(title, 
                  do.call(plot_grid, c(biplots.final, ncol = ncol, nrow = nrow)), 
                  ncol = 1, rel_heights = c(0.1, 1))
      }
    }

  # triangular layout?
  } else if (!triangle) {
    if (is.null(ncol)) {
      ncol <- length(components) - 1
    } else {
      ncol <- ncol
    }

    if (is.null(nrow)) {
      nrow <- ceiling(ncol / 2) + 1
    } else {
      nrow <- nrow
    }

    # remove axis labels and ticks?
    if (!plotaxes) {
      margin <- margin + theme(
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title = element_blank())
    } else if (plotaxes) {
      margin <- margin
    }

    # apply titles / axis changes
    biplots <- lapply(biplots, '+', margin)

    # return plot?
    if (returnPlot) {
      if(shared.legend) {
        legend <- get_legend(biplots[[2]] + theme(legend.box.margin = margin(0, 0, 0, 12)))
        bigrid <- plot_grid(title, 
                            do.call(plot_grid, c(lapply(biplots, "+", theme(legend.position = "none")), ncol = ncol, nrow = nrow)),
                            ncol = 1, rel_heights = c(0.1, 1))
        return(plot_grid(legend, bigrid, rel_widths = c(0.5, ncol)))
      } else {
        return(plot_grid(title, 
                         do.call(plot_grid, c(biplots, ncol = ncol, nrow = nrow)), 
                         ncol = 1, rel_heights = c(0.1, 1)))
      }
    } else if (!returnPlot) {
      if(shared.legend) {
        legend <- get_legend(biplots[[2]] + theme(legend.box.margin = margin(0, 0, 0, 12)))
        bigrid <- plot_grid(title, 
                            do.call(plot_grid, c(lapply(biplots, "+", theme(legend.position = "none")), ncol = ncol, nrow = nrow)), 
                            ncol = 1, rel_heights = c(0.1, 1))
        plot_grid(legend, bigrid, rel_widths = c(0.5, ncol))
      } else {
        plot_grid(title, 
                  do.call(plot_grid, c(biplots, ncol = ncol, nrow = nrow)), 
                  ncol = 1, rel_heights = c(0.1, 1))
      }
    }
  }
}
