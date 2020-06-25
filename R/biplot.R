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
  drawConnectorsLoadings = TRUE,
  widthConnectorsLoadings = 0.5,
  colConnectorsLoadings = 'grey50',
  lengthLoadingsArrowsFactor = 1.5,
  colLoadingsArrows = 'black',
  widthLoadingsArrows = 0.5,
  alphaLoadingsArrow = 1.0,
  colby = NULL,
  colkey = NULL,
  singlecol = NULL,
  shape = NULL,
  shapekey = NULL,
  pointSize = 3.0,
  legendPosition = 'none',
  legendLabSize = 12,
  legendIconSize = 5.0,
  encircleByGroup = FALSE,
  encircleFill = TRUE,
  encircleAlpha = 1/4,
  encircleLineSize = 0.25,
  ellipse = FALSE,
  ellipseConf = 0.95,
  ellipseFill = TRUE,
  ellipseAlpha = 1/4,
  ellipseLineSize = 0.25,
  xlim = c(min(pcaobj$rotated[,x]) - 2, max(pcaobj$rotated[,x]) + 2),
  ylim = c(min(pcaobj$rotated[,y]) - 2, max(pcaobj$rotated[,y]) + 2),
  lab = rownames(pcaobj$metadata),
  labSize = 3.0,
  labhjust = 1.5,
  labvjust = 0,
  boxedLabels = FALSE,
  selectLab = NULL,
  drawConnectors = TRUE,
  widthConnectors = 0.5,
  colConnectors = 'grey50',
  xlab = paste0(x, ', ',
    round(pcaobj$variance[x], digits = 2),
    '% variation'),
  xlabAngle = 0,
  xlabhjust = 0.5,
  xlabvjust = 0.5,
  ylab = paste0(y, ', ',
    round(pcaobj$variance[y], digits = 2),
    '% variation'),
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
      legend.title = element_blank())

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
      rownames(pcaobj$loadings)[xidx][seq_len(1:ntopLoadings)],
      rownames(pcaobj$loadings)[yidx][seq_len(1:ntopLoadings)]))

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
          plot <- plot + coord_equal() +
            geom_label_repel(data = pcaobj$loadings[vars,], 
              aes(label = vars,
                x = pcaobj$loadings[vars,x] * r * lengthLoadingsArrowsFactor,
                y = pcaobj$loadings[vars,y] * r * lengthLoadingsArrowsFactor,
              hjust = 0),
              color = colLoadingsNames,
              size = sizeLoadingsNames,
              segment.color = colConnectorsLoadings,
              segment.size = widthConnectorsLoadings)
        } else {
          plot <- plot + coord_equal() +
            geom_text_repel(data = pcaobj$loadings[vars,], 
              aes(label = vars,
                x = pcaobj$loadings[vars,x] * r * lengthLoadingsArrowsFactor,
                y = pcaobj$loadings[vars,y] * r * lengthLoadingsArrowsFactor,
              hjust = 0),
              color = colLoadingsNames,
              size = sizeLoadingsNames,
              segment.color = colConnectorsLoadings,
              segment.size = widthConnectorsLoadings)
        }
      } else {
        if (boxedLoadingsNames) {
          plot <- plot + coord_equal() +
            geom_label(data = pcaobj$loadings[vars,], 
              aes(label = vars,
                x = pcaobj$loadings[vars,x] * r * lengthLoadingsArrowsFactor,
                y = pcaobj$loadings[vars,y] * r * lengthLoadingsArrowsFactor,
              hjust = 0),
              color = colLoadingsNames,
              size = sizeLoadingsNames)
        } else {
          plot <- plot + coord_equal() +
            geom_text(data = pcaobj$loadings[vars,], 
              aes(label = vars,
                x = pcaobj$loadings[vars,x] * r * lengthLoadingsArrowsFactor,
                y = pcaobj$loadings[vars,y] * r * lengthLoadingsArrowsFactor,
              hjust = 0),
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
          size = labSize,
          segment.color = colConnectors,
          segment.size = widthConnectors,
          hjust = labhjust,
          vjust = labvjust)
    } else if (drawConnectors && !is.null(selectLab)) {
      plot <- plot + labFun(
        data=subset(plotobj,
          !is.na(plotobj[,'lab'])),
          aes(label = lab),
          size = labSize,
          segment.color = colConnectors,
          segment.size = widthConnectors,
          hjust = labhjust,
          vjust = labvjust)
    } else if (!drawConnectors && !is.null(selectLab)) {
      if (boxedLabels) {
        plot <- plot + labFun(
          data=subset(plotobj,
            !is.na(plotobj[,'lab'])),
            aes(label = lab),
            size = labSize,
            hjust = labhjust,
            vjust = labvjust)
      } else {
        plot <- plot + labFun(
          data=subset(plotobj,
            !is.na(plotobj[,'lab'])),
            aes(label = lab),
            size = labSize,
            check_overlap = TRUE,
            hjust = labhjust,
            vjust = labvjust)
      }
    } else if (!drawConnectors && is.null(selectLab)) {
      if (boxedLabels) {
        plot <- plot + labFun(
          data = plotobj,
            aes(label = lab),
            size = labSize,
            check_overlap = TRUE,
            hjust = labhjust,
            vjust = labvjust)
      } else {
        plot <- plot + labFun(
          data = plotobj,
            aes(label = lab),
            size = labSize,
            check_overlap = TRUE,
            hjust = labhjust,
            vjust = labvjust)
      }
    }
  }

  # encircle
  if (encircleByGroup) {
    if (encircleFill) {
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
            colour = col),
          alpha = encircleAlpha,
          size = encircleLineSize,
          show.legend = FALSE,
          na.rm = TRUE)
    }
  }

  # ellipse
  if (ellipse) {
    if (ellipseFill) {
      plot <- plot +
        stat_ellipse(
          aes(group = col,
            fill = col),
          colour = 'black',
          geom = 'polygon',
          level = ellipseConf,
          alpha = ellipseAlpha,
          size = ellipseLineSize,
          show.legend = FALSE,
          na.rm = TRUE)
    } else {
      plot <- plot +
        stat_ellipse(
          aes(group = col),
          colour = 'black',
          geom = 'polygon',
          level = ellipseConf,
          alpha = ellipseAlpha,
          size = ellipseLineSize,
          show.legend = FALSE,
          na.rm = TRUE)
    }
  }

  # return plot?
  if (returnPlot) {
    return(plot)
  } else if (!returnPlot) {
    plot
  }
}

