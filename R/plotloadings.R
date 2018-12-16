plotloadings <- function(
  pcaobj,
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
{
  # filter in the variables in the top percentRetain of the loadings range
  x <- pcaobj$loadings[,components]
  retain <- c()
  for (i in seq(length(components))) {
    # for each PC, based on the loadings range, calculate the rangeRetain
    # fraction of the range
    offset <- (max(x[,i]) - min(x[,i])) * rangeRetain

    # to obtain upper and lower cut-offs, offset max and min by the offset
    uppercutoff <- max(x[,i]) - offset
    lowercutoff <- min(x[,i]) + offset

    # build a consensus list of variables that will be included
    retain <- unique(c(retain,
      which((x[,i] >= uppercutoff) == TRUE),
      which((x[,i] <= lowercutoff) == TRUE)))
  }
  message('-- variables retained:')
  message(paste0(rownames(x)[retain], collapse = ', '))
  x <- x[retain,]

  PC <- Loading <- NULL

  # create a plot object (2-col df) of PC names and
  # explained variance of each
  x <- data.frame(rownames(x), x[,components])
  colnames(x)[1] <- 'var'
  plotobj <- melt(x, id = 'var')
  colnames(plotobj) <- c('var','PC','Loading')

  # convert loadings to absolute loadings?
  if (absolute == TRUE) {
    plotobj$Loading <- abs(plotobj$Loading)
  } else if (absolute == FALSE) {
    plotobj$Loading <- plotobj$Loading
  }

  # create a base theme that will later be modified
  th <- theme_bw(base_size=24) +

    theme(
      legend.background=element_rect(),

      plot.title=element_text(angle=0, size=titleLabSize,
        face='bold', vjust=1),
      plot.subtitle=element_text(angle = 0, size = subtitleLabSize,
        face = 'plain', vjust = 1),
      plot.caption=element_text(angle = 0, size = captionLabSize,
        face = 'plain', vjust = 1),

      axis.line = element_line(size=1.5, colour = 'black'),

      axis.text.x=element_text(angle = xlabAngle, size = axisLabSize,
        hjust = xlabhjust, vjust = xlabvjust),
      axis.text.y=element_text(angle = ylabAngle, size = axisLabSize,
        hjust = ylabhjust, vjust = ylabvjust),
      axis.title=element_text(size=axisLabSize),

      legend.position=legendPosition,
      legend.direction = 'horizontal',
      legend.box = 'horizontal',
      legend.key=element_blank(),
      legend.key.size=unit(0.5, 'cm'),
      legend.text=element_text(size=legendLabSize),

      title=element_text(size=legendLabSize),
      legend.title=element_blank())

  # create the plot object as geom_point
  plot <- ggplot(plotobj, aes(x = PC, y = Loading,
    size = Loading,
    fill = Loading)) + th +

    geom_point(shape = shape) +

    # add title, subtitle, caption
    labs(title = title, 
       subtitle = subtitle,
       caption = caption) +

    # add xy axis labeling
    labs(x = xlab, y = ylab, size = ylab, fill = ylab) +

    guides(fill = guide_legend(),
      size = guide_legend(),
      colour = guide_legend(override.aes = list(size = legendIconSize)))

  # scale the size of the geom_points based on the specified
  # shapeSizeRange, c(minSize, maxSize))
  plot <- plot + scale_size(range = shapeSizeRange)

  # colour the geom_points as a gradient based on the 2 or 3 colours
  # passed to 'col'
  if (length(col) == 2) {
    plot <- plot +
       scale_fill_continuous(low = col[1],
         high = col[2])
  } else if (length(col) == 3) {
    plot <- plot +
       scale_fill_gradient2(low = col[1],
         mid = col[2],
         high = col[3],
         midpoint = colMidpoint,
         space='Lab')
  }

  # add elements to the plot for xy axis limits
  if (!is.null(xlim)) {
    plot <- plot + xlim(xlim[1], xlim[2])
  }
  if (!is.null(ylim)) {
    plot <- plot + ylim(ylim[1], ylim[2])
  }

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

  # For labeling with geom_text_repel (connectors) and
  # geom_text(.., check_overlap = TRUE)
  if (drawConnectors == TRUE) {
    plot <- plot + geom_text_repel(
      data = plotobj,
        aes(label = as.character(var)),
        size = labSize,
        nudge_x = ifelse(positionConnectors == 'left', -0.75,
          ifelse(positionConnectors == 'right', 0.75, 0.0)),
        direction = 'y',
        segment.color = colConnectors,
        segment.size = widthConnectors,
        arrow = arrow(length = lengthConnectors,
          type = typeConnectors, ends = endsConnectors),
        show.legend = FALSE,
        hjust = labhjust,
        vjust = labvjust)
  } else if (drawConnectors == FALSE) {
    plot <- plot + geom_text(
      data = plotobj,
        aes(label = as.character(var)),
        size = labSize,
        check_overlap = TRUE,
        show.legend = FALSE,
        hjust = labhjust,
        vjust = labvjust)
  }

  # return plot?
  if (returnPlot == TRUE) {
    return(plot)
  } else if (returnPlot == FALSE) {
    plot
  }
}
