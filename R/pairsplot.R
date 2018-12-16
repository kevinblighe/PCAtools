pairsplot <- function(
  pcaobj,
  components = getComponents(pcaobj, seq(1, 5)),
  triangle = TRUE,
  trianglelabSize = 18,
  plottitles = FALSE,
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
  legendPosition = 'none',
  legendLabSize = 6,
  legendIconSize = 1.5,
  xlim = NULL,
  ylim = NULL,
  lab = FALSE,
  labSize = 1.5,
  labhjust = 1.5,
  labvjust = 0,
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
  title = '',
  subtitle = '',
  caption = '',
  titleLabSize = 10,
  subtitleLabSize = 6,
  captionLabSize = 6,
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

  # beginning of the master loop (contains nested loop)
  # biplots will be created on a pairwise basis
  for (i in seq_along(components)) {

    # if the triangular layout is specified, large titles are added
    # as ggdraw() objects
    if (triangle == TRUE) {
      nplots <- nplots + 1 # increment nplots

      biplots[[nplots]] <- ggdraw() +
        draw_label(x = 0.5, y = 0.5,
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
          title = title,
          subtitle = subtitle,
          caption = caption,
          titleLabSize = titleLabSize,
          subtitleLabSize = subtitleLabSize,
          captionLabSize = captionLabSize,
          pointSize = pointSize,
          labSize = labSize,
          labhjust = labhjust,
          labvjust = labvjust,
          legendPosition = legendPosition,
          legendLabSize = legendLabSize,
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

  # triangular layout?
  # with triangular layout, empty space is filled as empty
  # plots with label = ''
  if (triangle == TRUE) {
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

    # remove or keep titles, subtitles, caption
    # (removing saves valuable space)?
    if (plottitles == FALSE) {
      margin <- theme(
        plot.margin = margingaps,
        plot.title = element_blank(),
        plot.subtitle = element_blank(),
        plot.caption = element_blank())
    } else if (plottitles == TRUE) {
      margin <- theme(
        plot.margin = margingaps,
        plot.title = element_text(angle = 0, size = titleLabSize,
          face = 'bold', vjust = 1),
        plot.subtitle = element_text(angle = 0, size = subtitleLabSize,
          face = 'plain', vjust = 1),
        plot.caption=element_text(angle = 0, size = captionLabSize,
          face = 'plain', vjust = 1))
    }

    # remove axis labels and ticks?
    if (plotaxes == FALSE) {
      margin <- margin + theme(
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title = element_blank())
    } else if (plotaxes == TRUE) {
      margin <- margin
    }

    # apply titles / axis changes
    biplots.final <- lapply(biplots.final, '+', margin)
    biplots.final <- lapply(biplots.final, '+', coord_flip())

    # return plot?
    if (returnPlot == TRUE) {
      return(do.call(plot_grid, c(biplots.final, ncol = ncol, nrow = nrow)))
    } else if (returnPlot == FALSE) {
      do.call(plot_grid, c(biplots.final, ncol = ncol, nrow = nrow))
    }

  # triangular layout?
  } else if (triangle == FALSE) {
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

    # if triangle == FALSE, titles, subtitles, captions
    # are not included
    margin <- theme(
      plot.margin = margingaps,
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      plot.caption = element_blank())

    # remove axis labels and ticks?
    if (plotaxes == FALSE) {
      margin <- margin + theme(
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title = element_blank())
    } else if (plotaxes == TRUE) {
      margin <- margin
    }

    # apply titles / axis changes
    biplots <- lapply(biplots, '+', margin)

    # return plot?
    if (returnPlot == TRUE) {
      return(do.call(plot_grid, c(biplots, ncol = ncol, nrow = nrow)))
    } else if (returnPlot == FALSE) {
      do.call(plot_grid, c(biplots, ncol = ncol, nrow = nrow))
    }
  }
}
