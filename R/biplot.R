#' Draw a bi-plot, comparing 2 selected principal components / eigenvectors.
#'
#' @param pcaobj Object of class 'pca' created by pca().
#' @param x A principal component to plot on x-axis. All principal component
#'   names are stored in pcaobj$label.
#' @param y A principal component to plot on y-axis. All principal component
#'   names are stored in pcaobj$label.
#' @param showLoadings Logical, indicating whether or not to overlay
#'   variable loadings.
#' @param ntopLoadings If showLoadings == TRUE, select this many variables
#'   based on absolute ordered variable loading for each PC in the biplot.
#'   As a result of looking across 2 PCs, it can occur whereby greater than
#'   this number are actually displayed.
#' @param showLoadingsNames Logical, indicating to show variable loadings names
#'   or not.
#' @param colLoadingsNames If 'showLoadings == TRUE', colour of text labels.
#' @param sizeLoadingsNames If 'showLoadings == TRUE', size of text labels.
#' @param boxedLoadingsNames Logical, if 'showLoadings == TRUE', draw text
#'   labels in boxes.
#' @param fillBoxedLoadings When 'boxedLoadingsNames == TRUE', this controls
#'   the background fill of the boxes. To control both the fill and
#'   transparency, user can specify a value of the form
#'   'alpha(<colour>, <alpha>)'.
#' @param drawConnectorsLoadings If 'showLoadings == TRUE', draw line connectors
#'   to the variable loadings arrows in order to fit more labels in the plot
#'   space.
#' @param widthConnectorsLoadings If 'showLoadings == TRUE', width of the line
#'   connectors drawn to the variable loadings arrows.
#' @param colConnectorsLoadings If 'showLoadings == TRUE', colour of the line
#'   connectors drawn to the variable loadings arrows.
#' @param lengthLoadingsArrowsFactor If 'showLoadings == TRUE', multiply the
#'   internally-determined length of the variable loadings arrows by this
#'   factor.
#' @param colLoadingsArrows If showLoadings == TRUE, colour of the variable
#'   loadings arrows.
#' @param widthLoadingsArrows If showLoadings == TRUE, width of the variable
#'   loadings arrows.
#' @param alphaLoadingsArrow If showLoadings == TRUE, colour transparency of
#'   the variable loadings arrows.
#' @param colby If NULL, all points will be coloured differently. If not NULL,
#'   value is assumed to be a column name in pcaobj$metadata relating to some
#'   grouping/categorical variable.
#' @param colkey Vector of name-value pairs relating to value passed to 'col',
#'   e.g., c(A='forestgreen', B='gold').
#' @param colLegendTitle Title of the legend for the variable specified
#'   by 'colby'.
#' @param singlecol If specified, all points will be shaded by this colour.
#'   Overrides 'col'.
#' @param shape If NULL, all points will be have the same shape. If not NULL,
#'   value is assumed to be a column name in pcaobj$metadata relating to some
#'   grouping/categorical variable.
#' @param shapekey Vector of name-value pairs relating to value passed to
#'   'shape', e.g., c(A=10, B=21).
#' @param shapeLegendTitle Title of the legend for the variable specified
#'   by 'shape'.
#' @param pointSize Size of plotted points.
#' @param legendPosition Position of legend ('top', 'bottom', 'left', 'right',
#'   'none').
#' @param legendLabSize Size of plot legend text.
#' @param legendTitleSize Size of plot legend title text.
#' @param legendIconSize Size of plot legend icons / symbols.
#' @param encircle Logical, indicating whether to draw a polygon around
#'   the groups specified by 'colby'.
#' @param encircleFill Logical, if 'encircle == TRUE', this determines
#'   whether to fill the encircled region or not.
#' @param encircleFillKey Vector of name-value pairs relating to value passed to
#'   'encircleFill', e.g., c(A='forestgreen', B='gold'). If NULL, the fill
#'   is controlled by whatever has already been used for 'colby' / 'colkey'.
#' @param encircleAlpha Alpha for purposes of controlling colour transparency of
#'   the encircled region. Used when 'encircle == TRUE'.
#' @param encircleLineSize Line width of the encircled line when
#'   'encircle == TRUE'.
#' @param encircleLineCol Colour of the encircled line when
#'   'encircle == TRUE'.
#' @param ellipse Logical, indicating whether to draw a data ellipse around
#'   the groups specified by 'colby'.
#' @param ellipseType [paraphrased from
#'   https://ggplot2.tidyverse.org/reference/stat_ellipse.html]
#'   The type of ellipse. "t" assumes a multivariate t-distribution, while
#'   "norm" assumes a multivariate normal distribution. "euclid" draws a circle with
#'   the radius equal to level, representing the euclidean distance from the center.
#'   This ellipse probably won't appear circular unless coord_fixed() is applied.
#' @param ellipseLevel [paraphrased from
#'   https://ggplot2.tidyverse.org/reference/stat_ellipse.html]
#'   The level at which to draw an ellipse, or, if ellipseType="euclid", the radius of the circle to be drawn.
#' @param ellipseSegments [from
#'   https://ggplot2.tidyverse.org/reference/stat_ellipse.html]
#'   The number of segments to be used in drawing the ellipse.

#' @param ellipseFill Logical, if 'ellipse == TRUE', this determines
#'   whether to fill the region or not.
#' @param ellipseFillKey Vector of name-value pairs relating to value passed to
#'   'ellipseFill', e.g., c(A='forestgreen', B='gold'). If NULL, the fill
#'   is controlled by whatever has already been used for 'colby' / 'colkey'.
#' @param ellipseAlpha Alpha for purposes of controlling colour transparency of
#'   the ellipse region. Used when 'ellipse == TRUE'.
#' @param ellipseLineSize Line width of the ellipse line when 'ellipse == TRUE'.
#' @param ellipseLineCol Colour of the ellipse line when 'ellipse == TRUE'.
#' @param xlim Limits of the x-axis.
#' @param ylim Limits of the y-axis.
#' @param lab A vector containing labels to add to the plot. 
#' @param labSize Size of labels.
#' @param boxedLabels Logical, draw text labels in boxes.
#' @param selectLab A vector containing a subset of lab to plot.
#' @param drawConnectors Logical, indicating whether or not to connect plot
#'   labels to their corresponding points by line connectors.
#' @param widthConnectors Line width of connectors.
#' @param colConnectors Line colour of connectors.
#' @param maxoverlapsConnectors Equivalent of max.overlaps in ggrepel. Set to
#'   'Inf' to always display all labels when drawConnectors = TRUE.
#' @param directionConnectors direction in which to draw connectors.
#'   'both', 'x', or 'y'.
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
#' @details Draw a bi-plot, comparing 2 selected principal components / eigenvectors.
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
#'   rexp(col*row, rate = 0.1),
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
#'   biplot(p)
#'
#'   biplot(p, colby = 'Group', shape = 'Group')
#'
#'   biplot(p, colby = 'Group', colkey = c(A = 'forestgreen', B = 'gold'),
#'     legendPosition = 'right')
#'
#'   biplot(p, colby = 'Group', colkey = c(A='forestgreen', B='gold'),
#'     shape = 'Group', shapekey = c(A=10, B=21), legendPosition = 'bottom')
#'
#' @import ggplot2
#' @import ggrepel
#' 
#' @export
biplot <- function(
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
  ellipseType = 't',
  ellipseLevel = 0.95,
  ellipseSegments = 51,
  ellipseFill = TRUE,
  ellipseFillKey = NULL,
  ellipseAlpha = 1/4,
  ellipseLineSize = 0.25,
  ellipseLineCol = NULL,
  xlim = if(showLoadings || ellipse) c(
    min(pcaobj$rotated[,x]) - abs((min(pcaobj$rotated[,x])/100)*25),
    max(pcaobj$rotated[,x]) + abs((min(pcaobj$rotated[,x])/100)*25)) else c(
    min(pcaobj$rotated[,x]) - abs((min(pcaobj$rotated[,x])/100)*10),
    max(pcaobj$rotated[,x]) + abs((min(pcaobj$rotated[,x])/100)*10)),
  ylim = if(showLoadings || ellipse) c(
    min(pcaobj$rotated[,y]) - abs((min(pcaobj$rotated[,y])/100)*25),
    max(pcaobj$rotated[,y]) + abs((min(pcaobj$rotated[,y])/100)*25)) else c(
    min(pcaobj$rotated[,y]) - abs((min(pcaobj$rotated[,y])/100)*10),
    max(pcaobj$rotated[,y]) + abs((min(pcaobj$rotated[,y])/100)*10)),
  lab = rownames(pcaobj$metadata),
  labSize = 3.0,
  boxedLabels = FALSE,
  selectLab = NULL,
  drawConnectors = TRUE,
  widthConnectors = 0.5,
  colConnectors = 'grey50',
  maxoverlapsConnectors = 15,
  directionConnectors = 'both',
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
{

  labFun <- xidx <- yidx <- NULL

  # create a base theme that will later be modified
  th <- theme_bw(base_size = 24) +

    theme(
      legend.background = element_rect(),

      plot.title = element_text(angle = 0, size = titleLabSize,
        face = 'bold', vjust = 1),
      plot.subtitle = element_text(angle = 0, size = subtitleLabSize,
        face = 'plain', vjust = 1),
      plot.caption = element_text(angle = 0, size = captionLabSize,
        face = 'plain', vjust = 1),

      axis.text.x = element_text(angle = xlabAngle, size = axisLabSize,
        hjust = xlabhjust, vjust = xlabvjust),
      axis.text.y = element_text(angle = ylabAngle, size = axisLabSize,
        hjust = ylabhjust, vjust = ylabvjust),
      axis.title = element_text(size=axisLabSize),

      legend.position = legendPosition,
      legend.key = element_blank(),
      legend.key.size = unit(0.5, 'cm'),
      legend.text = element_text(size = legendLabSize),

      title = element_text(size = legendLabSize),
      legend.title = element_text(size = legendTitleSize))

  # set plot data labels (e.g. sample names)
  plotobj <- NULL
  plotobj$x <- pcaobj$rotated[,x]
  plotobj$y <- pcaobj$rotated[,y]
  if (!is.null(lab)) {
    plotobj$lab <- lab
  }
  plotobj <- as.data.frame(plotobj, stringsAsFactors = FALSE)

  # If user has supplied values in selectLab, convert labels to
  # NA and then re-set with those in selectLab
  if (!is.null(selectLab)) {
    if (is.null(lab)) {
      stop(paste0('You have specified lab as NULL ',
        '- no labels can be selected!'))
    } else {
      names.new <- rep(NA, length(plotobj$lab))
      indices <- which(plotobj$lab %in% selectLab)
      names.new[indices] <- plotobj$lab[indices]
      plotobj$lab <- names.new
    }
  }

  # decide on how to colour the points, and specify the shape of these
  if (is.null(colby)) {
    if (!is.null(lab)) {
      plotobj$col <- lab
    } else {
      plotobj$col <- seq_len(length(pcaobj$yvars))
    }
  } else {
    plotobj$col <- pcaobj$metadata[,colby]
  }
  if (!is.null(shape)) {
    plotobj$shape <- pcaobj$metadata[,shape]
  }

  # create the plot object
  plot <- ggplot(plotobj, aes(x = x, y = y)) + th +

    guides(fill = guide_legend(),
      shape = guide_legend(),
      colour = guide_legend(override.aes = list(size = legendIconSize)))

  # if user specified a colour with 'singlecol', colour all points by this
  # otherwise, colour all points differently using ggplot engine.
  # shape of points remains independent of colouring
  if (is.null(singlecol)) {
    if (!is.null(shape)) {
      plot <- plot + geom_point(aes(color = col, shape = shape),
        size = pointSize)
    } else {
      plot <- plot + geom_point(aes(color = col),
        size = pointSize)
    }
  } else if (!is.null(singlecol)) {
    if (!is.null(shape)) {
      plot <- plot + geom_point(aes(color = singlecol, shape = shape),
        size = pointSize)
    } else {
      plot <- plot + geom_point(aes(color = singlecol),
        size = pointSize)
    }
  }

  # sort out custom colour pairing, and custom shapes
  if (!is.null(colkey)) {
    plot <- plot + scale_colour_discrete('') +
      scale_color_manual(values = colkey)
  }
  if (!is.null(shapekey)) {
    plot <- plot + scale_shape_manual(values = shapekey)
  }

  # plot loadings arrows?
  if (showLoadings) {
    # get top ntopLoadings to display
    xidx <- order(abs(pcaobj$loadings[,x]), decreasing = TRUE)
    yidx <- order(abs(pcaobj$loadings[,y]), decreasing = TRUE)
    vars <- unique(c(
      rownames(pcaobj$loadings)[xidx][seq_len(ntopLoadings)],
      rownames(pcaobj$loadings)[yidx][seq_len(ntopLoadings)]))

    # get scaling parameter to match between variable loadings and rotated loadings
    r <- min(
      (max(pcaobj$rotated[,x]) - min(pcaobj$rotated[,x]) /
        (max(pcaobj$loadings[,x]) - min(pcaobj$loadings[,x]))),
      (max(pcaobj$rotated[,y]) - min(pcaobj$rotated[,y]) /
        (max(pcaobj$loadings[,y]) - min(pcaobj$loadings[,y]))))

    plot <- plot +
      geom_segment(data = pcaobj$loadings[vars,],
        aes(x = 0, y = 0,
          xend = pcaobj$loadings[vars,x] * r * lengthLoadingsArrowsFactor,
          yend = pcaobj$loadings[vars,y] * r * lengthLoadingsArrowsFactor),
        arrow = arrow(length = unit(1/2, 'picas'), ends = 'last'), 
        color = colLoadingsArrows,
        size = widthLoadingsArrows,
        alpha = alphaLoadingsArrow,
        show.legend = NA)

    if (showLoadingsNames) {
      if (drawConnectorsLoadings) {
        if (boxedLoadingsNames) {
          plot <- plot +
            geom_label_repel(data = pcaobj$loadings[vars,], 
              aes(label = vars,
                x = pcaobj$loadings[vars,x] * r * lengthLoadingsArrowsFactor,
                y = pcaobj$loadings[vars,y] * r * lengthLoadingsArrowsFactor),
              xlim = c(NA, NA),
              ylim = c(NA, NA),
              color = colLoadingsNames,
              size = sizeLoadingsNames,
              fill = fillBoxedLoadings,
              segment.color = colConnectorsLoadings,
              segment.size = widthConnectorsLoadings,
              direction = directionConnectors,
              max.overlaps = maxoverlapsConnectors)
        } else {
          plot <- plot +
            geom_text_repel(data = pcaobj$loadings[vars,], 
              aes(label = vars,
                x = pcaobj$loadings[vars,x] * r * lengthLoadingsArrowsFactor,
                y = pcaobj$loadings[vars,y] * r * lengthLoadingsArrowsFactor),
              xlim = c(NA, NA),
              ylim = c(NA, NA),
              color = colLoadingsNames,
              size = sizeLoadingsNames,
              segment.color = colConnectorsLoadings,
              segment.size = widthConnectorsLoadings,
              direction = directionConnectors,
              max.overlaps = maxoverlapsConnectors)
        }
      } else {
        if (boxedLoadingsNames) {
          plot <- plot +
            geom_label(data = pcaobj$loadings[vars,], 
              aes(label = vars,
                x = pcaobj$loadings[vars,x] * r * lengthLoadingsArrowsFactor,
                y = pcaobj$loadings[vars,y] * r * lengthLoadingsArrowsFactor),
              color = colLoadingsNames,
              size = sizeLoadingsNames,
              fill = NA)
        } else {
          plot <- plot +
            geom_text(data = pcaobj$loadings[vars,], 
              aes(label = vars,
                x = pcaobj$loadings[vars,x] * r * lengthLoadingsArrowsFactor,
                y = pcaobj$loadings[vars,y] * r * lengthLoadingsArrowsFactor),
              color = colLoadingsNames,
              size = sizeLoadingsNames,
              check_overlap = TRUE)
        }
      }
    }
  }

  # add elements to the plot for xy labeling and axis limits
  plot <- plot + xlab(xlab) + ylab(ylab)
  if (!is.null(xlim)) {
    plot <- plot + xlim(xlim[1], xlim[2])
  }
  if (!is.null(ylim)) {
    plot <- plot + ylim(ylim[1], ylim[2])
  }

  # add elements to the plot for title, subtitle, caption, and legend titles
  plot <- plot + labs(title = title, 
    subtitle = subtitle, caption = caption,
    fill = '', colour = colLegendTitle, shape = shapeLegendTitle)

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

  # labeling
  if (boxedLabels) {
    if (drawConnectors) {
      labFun <- function(...) geom_label_repel(...)
    } else {
      labFun <- function(...) geom_label(...)
    }
  } else {
    if (drawConnectors) {
      labFun <- function(...) geom_text_repel(...)
    } else {
      labFun <- function(...) geom_text(...)
    }
  }

  if (!is.null(lab)) {
    if (drawConnectors && is.null(selectLab)) {
      plot <- plot + labFun(
        data = plotobj,
          aes(label = lab),
          xlim = c(NA, NA),
          ylim = c(NA, NA),
          size = labSize,
          segment.color = colConnectors,
          segment.size = widthConnectors,
          direction = directionConnectors,
          max.overlaps = maxoverlapsConnectors)
    } else if (drawConnectors && !is.null(selectLab)) {
      plot <- plot + labFun(
        data=subset(plotobj,
          !is.na(plotobj[,'lab'])),
          aes(label = lab),
          xlim = c(NA, NA),
          ylim = c(NA, NA),
          size = labSize,
          segment.color = colConnectors,
          segment.size = widthConnectors,
          direction = directionConnectors,
          max.overlaps = maxoverlapsConnectors)
    } else if (!drawConnectors && !is.null(selectLab)) {
      if (boxedLabels) {
        plot <- plot + labFun(
          data=subset(plotobj,
            !is.na(plotobj[,'lab'])),
            aes(label = lab),
            xlim = c(NA, NA),
            ylim = c(NA, NA),
            size = labSize)
      } else {
        plot <- plot + labFun(
          data=subset(plotobj,
            !is.na(plotobj[,'lab'])),
            aes(label = lab),
            xlim = c(NA, NA),
            ylim = c(NA, NA),
            size = labSize,
            check_overlap = TRUE)
      }
    } else if (!drawConnectors && is.null(selectLab)) {
      if (boxedLabels) {
        plot <- plot + labFun(
          data = plotobj,
            aes(label = lab),
            xlim = c(NA, NA),
            ylim = c(NA, NA),
            size = labSize,
            check_overlap = TRUE)
      } else {
        plot <- plot + labFun(
          data = plotobj,
            aes(label = lab),
            xlim = c(NA, NA),
            ylim = c(NA, NA),
            size = labSize,
            check_overlap = TRUE)
      }
    }
  }

  # encircle
  if (encircle) {
    if (encircleFill) {
      if (is.null(encircleLineCol)) {
        plot <- plot +
          ggalt::geom_encircle(
            aes(group = col,
              fill = col,
              colour = col),
            alpha = encircleAlpha,
            size = encircleLineSize,
            show.legend = FALSE,
            na.rm = TRUE)
      } else {
        plot <- plot +
          ggalt::geom_encircle(
            aes(group = col,
              fill = col),
            colour = encircleLineCol,
            alpha = encircleAlpha,
            size = encircleLineSize,
            show.legend = FALSE,
            na.rm = TRUE)
      }
    } else {
      if (is.null(encircleLineCol)) {
        plot <- plot +
          ggalt::geom_encircle(
            aes(group = col,
              colour = col),
            fill = NA,
            alpha = encircleAlpha,
            size = encircleLineSize,
            show.legend = FALSE,
            na.rm = TRUE)
      } else {
        plot <- plot +
          ggalt::geom_encircle(
            aes(group = col),
            colour = encircleLineCol,
            fill = NA,
            alpha = encircleAlpha,
            size = encircleLineSize,
            show.legend = FALSE,
            na.rm = TRUE)
      }
    }

    if (encircleFill) {
      if (is.null(encircleFillKey)) {
        if (!is.null(colkey)) {
          plot <- plot + scale_fill_manual(values = colkey)
        }
      } else {
          plot <- plot + scale_fill_manual(values = encircleFillKey)
      }
    }
  }

  # ellipse
  if (ellipse) {
    if (ellipseFill) {
      if (is.null(ellipseLineCol)) {
        plot <- plot +
          stat_ellipse(
            aes(group = col,
              fill = col,
              colour = col),
            geom = 'polygon',
            type = ellipseType,
            level = ellipseLevel,
            segments = ellipseSegments,
            alpha = ellipseAlpha,
            size = ellipseLineSize,
            show.legend = FALSE,
            na.rm = TRUE)
      } else {
        plot <- plot +
          stat_ellipse(
            aes(group = col,
              fill = col),
            colour = ellipseLineCol,
            geom = 'polygon',
            type = ellipseType,
            level = ellipseLevel,
            segments = ellipseSegments,
            alpha = ellipseAlpha,
            size = ellipseLineSize,
            show.legend = FALSE,
            na.rm = TRUE)
      }
    } else {
      if (is.null(ellipseLineCol)) {
        plot <- plot +
          stat_ellipse(
            aes(group = col,
              colour = col),
            fill = NA,
            geom = 'polygon',
            type = ellipseType,
            level = ellipseLevel,
            segments = ellipseSegments,
            alpha = ellipseAlpha,
            size = ellipseLineSize,
            show.legend = FALSE,
            na.rm = TRUE)
      } else {
        plot <- plot +
          stat_ellipse(
            aes(group = col),
            colour = ellipseLineCol,
            fill = NA,
            geom = 'polygon',
            type = ellipseType,
            level = ellipseLevel,
            segments = ellipseSegments,
            alpha = ellipseAlpha,
            size = ellipseLineSize,
            show.legend = FALSE,
            na.rm = TRUE)
      }
    }

    if (ellipseFill) {
      if (is.null(ellipseFillKey)) {
        if (!is.null(colkey)) {
          plot <- plot + scale_fill_manual(values = colkey)
        }
      } else {
          plot <- plot + scale_fill_manual(values = ellipseFillKey)
      }
    }
  }

  plot <- plot + coord_cartesian(clip = 'off')

  # return plot?
  if (returnPlot) {
    return(plot)
  } else if (!returnPlot) {
    plot
  }
}

