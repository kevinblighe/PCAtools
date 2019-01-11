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
