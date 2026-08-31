#' Draw a bi-plot, comparing 2 selected principal components / eigenvectors.
#'
#' @param pcaobj Object of class 'pca' created by [pca()].
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
#' @param LoadingsNames Names of variable loadings to display. If provided,
#'   this overrides 'ntopLoadings'.
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
#' @param max.overlaps Equivalent of max.overlaps in ggrepel. Set to
#'   'Inf' to always display all labels when drawConnectors = TRUE.
#' @param maxoverlapsConnectors See max.overlaps.
#' @param min.segment.length When drawConnectors = TRUE, specifies the minimum
#'   length of the connector line segments.
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
#' @param flip_axes Logical, indicating whether to swap the x- and y-axes.
#'
#' @details Draw a bi-plot, comparing 2 selected principal components / eigenvectors.
#'
#' @return A `ggplot2` object.
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
  LoadingsNames = NULL,
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
    min(pcaobj$rotated[,x]) - abs((min(pcaobj$rotated[,x])/100)*35),
    max(pcaobj$rotated[,x]) + abs((min(pcaobj$rotated[,x])/100)*35)) else c(
    min(pcaobj$rotated[,x]) - abs((min(pcaobj$rotated[,x])/100)*10),
    max(pcaobj$rotated[,x]) + abs((min(pcaobj$rotated[,x])/100)*10)),
  ylim = if(showLoadings || ellipse) c(
    min(pcaobj$rotated[,y]) - abs((min(pcaobj$rotated[,y])/100)*35),
    max(pcaobj$rotated[,y]) + abs((min(pcaobj$rotated[,y])/100)*35)) else c(
    min(pcaobj$rotated[,y]) - abs((min(pcaobj$rotated[,y])/100)*10),
    max(pcaobj$rotated[,y]) + abs((min(pcaobj$rotated[,y])/100)*10)),
  lab = rownames(pcaobj$metadata),
  labSize = 3.0,
  boxedLabels = FALSE,
  selectLab = NULL,
  drawConnectors = TRUE,
  widthConnectors = 0.5,
  colConnectors = 'grey50',
  max.overlaps = 15,
  maxoverlapsConnectors = NULL,
  min.segment.length = 0,
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
  returnPlot = TRUE,
  flip_axes = FALSE)
{

  labFun <- xidx <- yidx <- NULL

  if (!is.null(maxoverlapsConnectors)) {
    max.overlaps <- maxoverlapsConnectors
  }

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
  xy_aes <- if (flip_axes) aes(x = y, y = x) else aes(x = x, y = y)
  plot <- ggplot(plotobj, xy_aes) + th +
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
    # get top ntopLoadings to display, or use user-supplied LoadingsNames
    xidx <- order(abs(pcaobj$loadings[,x]), decreasing = TRUE)
    yidx <- order(abs(pcaobj$loadings[,y]), decreasing = TRUE)
    if (!is.null(LoadingsNames)) {
      # validate provided names exist in loadings
      missing.names <- setdiff(LoadingsNames, rownames(pcaobj$loadings))
      if (length(missing.names) > 0) {
        stop(paste0('The following LoadingsNames are not present in pcaobj$loadings: ',
          paste(missing.names, collapse = ', ')))
      }
      vars <- unique(LoadingsNames)
    } else {
      vars <- unique(c(
        rownames(pcaobj$loadings)[xidx][seq_len(ntopLoadings)],
        rownames(pcaobj$loadings)[yidx][seq_len(ntopLoadings)]))
    }

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
              max.overlaps = max.overlaps,
              min.segment.length = min.segment.length)
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
              max.overlaps = max.overlaps,
              min.segment.length = min.segment.length)
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
    fill = NULL, colour = colLegendTitle, shape = shapeLegendTitle)

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
      linewidth = borderWidth))

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
          max.overlaps = max.overlaps,
          min.segment.length = min.segment.length)
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
          max.overlaps = max.overlaps,
          min.segment.length = min.segment.length)
    } else if (!drawConnectors && !is.null(selectLab)) {
      if (boxedLabels) {
        plot <- plot + labFun(
          data=subset(plotobj,
            !is.na(plotobj[,'lab'])),
            aes(label = lab),
            size = labSize)
      } else {
        plot <- plot + labFun(
          data=subset(plotobj,
            !is.na(plotobj[,'lab'])),
            aes(label = lab),

            size = labSize,
            check_overlap = TRUE)
      }
    } else if (!drawConnectors && is.null(selectLab)) {
      if (boxedLabels) {
        plot <- plot + labFun(
          data = plotobj,
            aes(label = lab),
            size = labSize,
            check_overlap = TRUE)
      } else {
        plot <- plot + labFun(
          data = plotobj,
            aes(label = lab),
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
          geom_encircle(
            aes(group = col,
              fill = col,
              colour = col),
            alpha = encircleAlpha,
            size = encircleLineSize,
            show.legend = FALSE,
            na.rm = TRUE)
        
      } else {
        plot <- plot +
          geom_encircle(
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
          geom_encircle(
            aes(group = col,
              colour = col),
            fill = NA,
            alpha = encircleAlpha,
            size = encircleLineSize,
            show.legend = FALSE,
            na.rm = TRUE)
        
        # plot <- plot +
        #     ggforce::geom_mark_ellipse(
        #         aes(group = col,
        #             colour = col),
        #         fill = NA,
        #         alpha = encircleAlpha,
        #         size = encircleLineSize,
        #         show.legend = FALSE,
        #         na.rm = TRUE)
      } else {
        plot <- plot +
          geom_encircle(
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

#' Custom Geom for Encircling Points in ggplot2
#' @importFrom ggplot2 ggproto Geom aes layer
#' @importFrom grid unit convertUnit get.gpar gpar xsplineGrob grobTree rectGrob
#' @importFrom grDevices chull
#' @importFrom scales alpha
#' @format NULL
#' @usage NULL
#' @author Jared Andrews, heavily based on ggalt code from Ben Bolker
#'   (\url{https://github.com/hrbrmstr/ggalt/blob/master/R/geom_encircle.r})
#' @export
GeomEncircle <- ggproto("GeomEncircle", Geom,
                        required_aes = c("x", "y"),
                        default_aes = aes(
                            colour = "black",
                            fill = NA,
                            alpha = 1,
                            linetype = 1,
                            size = 1,
                            s_shape = 0.5,
                            s_open = FALSE,
                            expand = 0.05,
                            spread = 0.1
                        ),
                        
    draw_group = function(data, panel_scales, coord) {
        # Apply coordinate transformation
        transformed_coords <- coord$transform(data, panel_scales)
        reference_row <- transformed_coords[1, , drop = FALSE]
        rownames(reference_row) <- NULL
        
        # Find center point of all coordinates
        center_point <- lapply(transformed_coords[, c("x", "y")], mean, na.rm = TRUE)
        
        # Get convex hull vertex indices
        hull_vertices <- grDevices::chull(transformed_coords[c("x", "y")])
        
        # Factory function for coordinate data frames
        build_coord_frame <- function(x_coords, y_coords) {
            non_xy_cols <- reference_row[!names(reference_row) %in% c("x", "y")]
            data.frame(x = x_coords, y = y_coords, non_xy_cols)
        }
        
        transformed_coords <- transformed_coords[hull_vertices, ]
        
        # Coordinate conversion utilities
        native_to_millimeters <- function(value, dimension = "x") {
            grid::convertUnit(
                grid::unit(value, "native"), "mm",
                typeFrom = "dimension", axisFrom = dimension, valueOnly = TRUE
            )
        }
        
        location_to_native <- function(value, dimension = "x") {
            grid::convertUnit(
                value, "native",
                typeFrom = "location", axisFrom = dimension, valueOnly = TRUE
            )
        }
        
        merge_native_snpc <- function(native_component, snpc_component, dimension = "x") {
            location_to_native(
                unit(native_component, "native") + unit(snpc_component, "snpc"),
                dir = dimension
            )
        }
        
        # Calculate unit vector from point2 to point1
        unit_vector <- function(point1, point2) {
            x_diff <- native_to_millimeters(point1$x - point2$x)
            y_diff <- native_to_millimeters(point1$y - point2$y)
            vector_length <- sqrt(x_diff^2 + y_diff^2)
            list(x = x_diff / vector_length, y = y_diff / vector_length)
        }
        
        # Handle edge cases for small point sets
        if (nrow(transformed_coords) == 1) {
            # Single point: expand into diamond
            transformed_coords <- with(transformed_coords, build_coord_frame(
                c(x, x + spread, x, x - spread),
                c(y + spread, y, y - spread, y)
            ))
        } else if (nrow(transformed_coords) == 2) {
            # Two points: create perpendicular diamond
            perpendicular_rotation <- matrix(c(0, 1, -1, 0), nrow = 2)
            direction <- unit_vector(transformed_coords[1, ], transformed_coords[2, ])
            offset <- c(perpendicular_rotation %*% unlist(direction)) * 
                transformed_coords$spread
            
            transformed_coords <- with(transformed_coords, {
                new_x <- c(x[1], center_point$x + offset[1], x[2], center_point$x - offset[1])
                new_y <- c(y[1], center_point$y + offset[2], y[2], center_point$y - offset[2])
                build_coord_frame(new_x, new_y)
            })
        }
        
        # Calculate outward directions from center
        outward_vectors <- unit_vector(transformed_coords, center_point)
        
        # Configure graphics parameters
        gpar_config <- grid::get.gpar()
        aesthetic_inputs <- c("colour", "linetype", "alpha", "fill", "size")
        gpar_outputs <- c("col", "lty", "alpha", "fill", "lwd")
        gpar_config[gpar_outputs] <- reference_row[aesthetic_inputs]
        
        # Generate the encircling spline
        grid::xsplineGrob(
            with(transformed_coords, 
                 unit(x, "npc") + outward_vectors$x * unit(expand, "snpc")),
            with(transformed_coords, 
                 unit(y, "npc") + outward_vectors$y * unit(expand, "snpc")),
            shape = transformed_coords$s_shape - 1,
            open = reference_row$s_open,
            gp = gpar_config
        )
    }
)

#' Automatically enclose points in a polygon
#'
#' Creates a smooth encircling polygon around a set of points using convex hull
#' calculation and xspline smoothing. Useful for highlighting groups of points
#' in scatter plots.
#'
#' @param mapping Set of aesthetic mappings created by \code{\link[ggplot2]{aes}}. 
#'   If specified and \code{inherit.aes = TRUE} (the default), it is combined with 
#'   the default mapping at the top level of the plot.
#' @param data The data to be displayed in this layer. If \code{NULL}, the default,
#'   the data is inherited from the plot data as specified in the call to 
#'   \code{\link[ggplot2]{ggplot}}.
#' @param stat The statistical transformation to use on the data for this layer,
#'   as a string.
#' @param position Position adjustment, either as a string, or the result of a call
#'   to a position adjustment function.
#' @param na.rm If \code{FALSE}, the default, missing values are removed with a warning.
#'   If \code{TRUE}, missing values are silently removed.
#' @param show.legend Logical. Should this layer be included in the legends?
#'   \code{NA}, the default, includes if any aesthetics are mapped.
#' @param inherit.aes If \code{FALSE}, overrides the default aesthetics, rather
#'   than combining with them.
#' @param ... Other arguments passed on to \code{\link[ggplot2]{layer}}. These are
#'   often aesthetics, used to set an aesthetic to a fixed value, like 
#'   \code{colour = "red"} or \code{size = 3}. They may also be parameters to the
#'   paired geom/stat. Additional parameters include:
#'   \describe{
#'     \item{s_shape}{Controls the shape of the spline (default = 0.5).}
#'     \item{s_open}{Logical indicating whether the spline should be open (default = FALSE).}
#'     \item{expand}{Amount to expand the encircling polygon outward (default = 0.05).}
#'     \item{spread}{Spread factor for single or double point sets (default = 0.1).}
#'   }
#'
#' @return A ggplot2 layer that can be added to a plot.
#'
#' @author Jared Andrews, heavily based on ggalt code from Ben Bolker
#'   (\url{https://github.com/hrbrmstr/ggalt/blob/master/R/geom_encircle.r})
#'
#' @export
#' @examples
#' \dontrun{
#' library(ggplot2)
#'
#' d <- data.frame(x=c(1,1,2),y=c(1,2,2)*100)
#'
#' gg <- ggplot(d,aes(x,y))
#' gg <- gg + scale_x_continuous(expand=c(0.5,1))
#' gg <- gg + scale_y_continuous(expand=c(0.5,1))
#'
#' gg + geom_encircle(s_shape=1, expand=0) + geom_point()
#'
#' gg + geom_encircle(s_shape=1, expand=0.1, colour="red") + geom_point()
#'
#' gg + geom_encircle(s_shape=0.5, expand=0.1, colour="purple") + geom_point()
#'
#' gg + geom_encircle(data=subset(d, x==1), colour="blue", spread=0.02) +
#'   geom_point()
#'
#' gg + geom_encircle(data=subset(d, x==2), colour="cyan", spread=0.04) +
#'   geom_point()
#'
#' gg <- ggplot(mpg, aes(displ, hwy))
#' gg + geom_encircle(data=subset(mpg, hwy>40)) + geom_point()
#' gg + geom_encircle(aes(group=manufacturer)) + geom_point()
#' gg + geom_encircle(aes(group=manufacturer,fill=manufacturer),alpha=0.4)+
#'        geom_point()
#' gg + geom_encircle(aes(group=manufacturer,colour=manufacturer))+
#'        geom_point()
#'
#' ss <- subset(mpg,hwy>31 & displ<2)
#'
#' gg + geom_encircle(data=ss, colour="blue", s_shape=0.9, expand=0.07) +
#'   geom_point() + geom_point(data=ss, colour="blue")
#' }
geom_encircle <- function(mapping = NULL, data = NULL, stat = "identity",
                          position = "identity", na.rm = FALSE, show.legend = NA,
                          inherit.aes = TRUE, ...) {
    layer(
        geom = GeomEncircle, mapping = mapping, data = data, stat = stat,
        position = position, show.legend = show.legend, inherit.aes = inherit.aes,
        params = list(na.rm = na.rm, ...)
    )
}
