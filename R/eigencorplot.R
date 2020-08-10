#' Correlate principal components to continuous variable metadata and test significancies of these.
#'
#' @param pcaobj Object of class 'pca' created by pca().
#' @param components The principal components to be included in the plot.
#' @param metavars A vector of column names in metadata representing continuos
#'   variables.
#' @param titleX X-axis title.
#' @param cexTitleX X-axis title cex.
#' @param rotTitleX X-axis title rotation in degrees.
#' @param colTitleX X-axis title colour.
#' @param fontTitleX X-axis title font style. 1, plain; 2, bold; 3, italic; 4,
#'   bold-italic.
#' @param titleY Y-axis title.
#' @param cexTitleY Y-axis title cex.
#' @param rotTitleY Y-axis title rotation in degrees.
#' @param colTitleY Y-axis title colour.
#' @param fontTitleY Y-axis title font style. 1, plain; 2, bold; 3, italic; 4,
#'   bold-italic.
#' @param cexLabX X-axis labels cex.
#' @param rotLabX X-axis labels rotation in degrees.
#' @param colLabX X-axis labels colour.
#' @param fontLabX X-axis labels font style. 1, plain; 2, bold; 3, italic; 4,
#'   bold-italic.
#' @param cexLabY Y-axis labels cex.
#' @param rotLabY Y-axis labels rotation in degrees.
#' @param colLabY Y-axis labels colour.
#' @param fontLabY Y-axis labels font style. 1, plain; 2, bold; 3, italic; 4,
#'   bold-italic.
#' @param posLab Positioning of the X- and Y-axis labels. 'bottomleft', bottom
#'   and left; 'topright', top and right; 'all', bottom / top and left /right;
#'   'none', no labels.
#' @param col Colour shade gradient for RColorBrewer.
#' @param posColKey Position of colour key. 'bottom', 'left', 'top', 'right'.
#' @param cexLabColKey Colour key labels cex.
#' @param cexCorval Correlation values cex.
#' @param colCorval Correlation values colour.
#' @param fontCorval Correlation values font style. 1, plain; 2, bold; 3,
#'   italic; 4, bold-italic.
#' @param scale Logical, indicating whether or not to scale the colour range
#'   to max and min cor values.
#' @param main Plot title.
#' @param cexMain Plot title cex.
#' @param rotMain Plot title rotation in degrees.
#' @param colMain Plot title colour.
#' @param fontMain Plot title font style. 1, plain; 2, bold; 3, italic; 4,
#'   bold-italic.
#' @param corFUN Correlation method: 'pearson', 'spearman', or 'kendall'.
#' @param corUSE Method for handling missing values (see documentation for cor
#'   function via ?cor). 'everything', 'all.obs', 'complete.obs',
#'   'na.or.complete', or 'pairwise.complete.obs'.
#' @param corMultipleTestCorrection Multiple testing p-value adjustment method.
#'   Any method from stats::p.adjust() can be used. Activating this function means that
#'   signifSymbols and signifCutpoints then relate to adjusted (not nominal) 
#'   p-values.
#' @param signifSymbols Statistical significance symbols to display beside
#'   correlation values.
#' @param signifCutpoints Cut-points for statistical significance.
#' @param colFrame Frame colour.
#' @param plotRsquared Logical, indicating whether or not to plot R-squared
#'   values.
#' @param returnPlot Logical, indicating whether or not to return the plot
#'   object.
#'
#' @details Correlate principal components to continuous variable metadata and test significancies of these.
#'
#' @return A \code{\link{lattice}} object.
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
#'   eigencorplot(p, components = getComponents(p, 1:10),
#'     metavars = c('ESR', 'CRP'))
#'
#' @importFrom stats cor cor.test p.adjust symnum
#' @importFrom lattice levelplot
#' @importFrom lattice panel.levelplot
#' @importFrom lattice ltext
#' @importFrom grDevices colorRampPalette
#' @importFrom methods is
#' 
#' @export
eigencorplot <- function(
  pcaobj,
  components = getComponents(pcaobj, seq_len(10)),
  metavars,
  titleX = '',
  cexTitleX = 1.0,
  rotTitleX = 0,
  colTitleX = 'black',
  fontTitleX = 2,
  titleY = '',
  cexTitleY = 1.0,
  rotTitleY = 0,
  colTitleY = 'black',
  fontTitleY = 2,
  cexLabX = 1.0,
  rotLabX = 0,
  colLabX = 'black',
  fontLabX = 2,
  cexLabY = 1.0,
  rotLabY = 0,
  colLabY = 'black',
  fontLabY = 2,
  posLab = 'bottomleft',
  col = c('blue4', 'blue3', 'blue2', 'blue1', 'white',
    'red1', 'red2', 'red3', 'red4'),
  posColKey = 'right',
  cexLabColKey = 1.0,
  cexCorval = 1.0,
  colCorval = 'black',
  fontCorval = 1,
  scale = TRUE,
  main = '',
  cexMain = 2,
  rotMain = 0,
  colMain = 'black',
  fontMain = 2,
  corFUN = 'pearson',
  corUSE = 'pairwise.complete.obs',
  corMultipleTestCorrection = 'none',
  signifSymbols = c('***', '**', '*', ''),
  signifCutpoints = c(0, 0.001, 0.01, 0.05, 1),
  colFrame = 'white',
  plotRsquared = FALSE,
  returnPlot = TRUE)

{
  data <- pcaobj$rotated
  metadata <- pcaobj$metadata

  # issue warning if any columns to use are not numeric --- This is kind of annoying
  # --- and there should be an option to force numeric conversion as anyway I'm doing it in advance
  
  for (i in seq_len(length(components))) {
    if(!is.numeric(data[,components[i]])) {
      warning(components[i],
        ' is not numeric - please check the source data',
        ' as everything will be converted to a matrix')
    }
  }
  for (i in seq_len(length(metavars))) {
    if(!is.numeric(metadata[,metavars[i]])) {
      warning(metavars[i],
        ' is not numeric - please check the source data',
        ' as non-numeric variables will be coerced to numeric')
    }
  }

  # convert the data for x and y to data matrix
  #	NAs are left NA
  #	Character (A-Z a-z) are converted to NA
  #	Character numbers are converted to integers
  #	Factors are converted to numbers based on level ordering
  xvals <- data.matrix(data[,which(colnames(data) %in% components)])
  yvals <- metadata[,which(colnames(metadata) %in% metavars)]

  ### code courtesy of aleighbrown
    # let's make sure that anything that isnt' numeric becomes a numeric
    # find all the character columns
    # ---select all the columns which are characters
    # using base R now for dependencies sake 
      chararcter_columns = unlist(lapply(yvals, is.numeric))  
      # negate it - basically if it
      chararcter_columns = !chararcter_columns
      # select only the names that are true 
      chararcter_columns = names(which(chararcter_columns))
      for (c in chararcter_columns) {
        yvals[, eval(quote(c))] = as.numeric(as.factor(yvals[, eval(quote(c))]))
      }
  ### END
    
  yvals <- data.matrix(yvals)

  # create correlation table
  corvals <- cor(xvals, yvals, use = corUSE, method = corFUN)

  # create a new df with same dimensions as corvals and fill with P values
  # total number of tests we perform
  N <- ncol(xvals) * ncol(yvals)
  pvals <- data.frame(pval = numeric(N),
                       i = numeric(N),
                       j = numeric(N))
  k <- 0
  for (i in seq_len(ncol(xvals))) {
    for (j in seq_len(ncol(yvals))) { 
      k <- k + 1
      pvals[k,'pval'] <- cor.test(xvals[,i],
                             yvals[,j],
                             use = corUSE,
                             method = corFUN)$p.value
      pvals[k,"i"] <- colnames(xvals)[i]
      pvals[k,"j"] <- colnames(yvals)[j]

    }
  }

  ### code courtesy of aleighbrown
    # -----if you want to adjust the p-values for multiple testing
    if(corMultipleTestCorrection != "none"){
      pvals$pval <- p.adjust(pvals$pval, method = corMultipleTestCorrection)
    }

    pvals <- reshape2::dcast(pvals, i ~ j, value.var = "pval")
    # -----make sure the pvals matchs the order of corrvals table
    rownames(pvals) <- pvals$i
    pvals$i <- NULL
    pvals <- pvals[match(rownames(corvals), rownames(pvals)), ]
    # ---make sure the columns are in the correct order
    pvals <- pvals[colnames(corvals)]
    # ------
  ### END

  # are we plotting R^2 values?
  if (plotRsquared==TRUE) {
    corvals <- corvals ^ 2
  }

  # determine max and min correlation values in order to define the range
  if (scale == FALSE && plotRsquared == TRUE) {
    iUpperRange <- 1
    iLowerRange <- 0
  } else if (scale == FALSE && plotRsquared == FALSE) {
    iUpperRange <- 1
    iLowerRange <- -1
  } else if (scale == TRUE) {
    max <- max(corvals)
    min <- min(corvals)
    if(abs(max) > abs(min)) {
      iUpperRange <- max + 0.01
      iLowerRange <- (max * (-1)) - 0.01
    } else {
      iUpperRange <- abs(min) + 0.01
      iLowerRange <- min - 0.01
    }
    if (plotRsquared==TRUE) {
      iUpperRange <- max + 0.1
      iLowerRange <- 0
    }
  }

  # define the colour scheme/palette
  cols <- colorRampPalette(col)

  # create a new df with same dimensions as corvals
  # fill with significances encoded with asterisks
  signif <- corvals
  for (i in seq_len(ncol(pvals))) {
    signif[,i] <- c(symnum(pvals[,i],
      corr = FALSE,
      na = FALSE,
      cutpoints = signifCutpoints,
      symbols = signifSymbols))
  }

  # create a new df with same dimensions as corvals
  # fill with r values merged with the encoded significances
  plotLabels <- corvals
  for (i in seq_len(nrow(corvals))) {
    for(j in seq_len(ncol(corvals))) {
      plotLabels[i,j] <- paste(round(corvals[i,j], 2),
        signif[i,j],
        sep='')
      colnames(plotLabels)[j] <- colnames(corvals)[j]
    }

    rownames(plotLabels)[i] <- rownames(corvals)[i]
  }

  # position of axis ticks
  if (posLab == 'bottomleft') {
    posLab = 1
    axisTicks = c(1,0)
  } else if (posLab == 'topright') {
    posLab = 2
    axisTicks = c(0,1)
  } else if (posLab == 'all') {
    posLab = 3
    axisTicks = c(1,1)
  } else if (posLab == 'none') {
    posLab = 0
    axisTicks = c(0,0)
  }

  # define a panel function for adding labels
  # labels are passed with z as a third dimension
  labels <- function(x, y, z, ...) {
    panel.levelplot(x, y, z, ...)
    ltext(x, y,
      labels = plotLabels,
      cex = cexCorval,
      col = colCorval,
      font = fontCorval)
  }

  # produce the levelplot
  l <- levelplot(
    data.matrix(corvals),
    xlab = list(label = titleX,
      cex = cexTitleX,
      rot = rotTitleX,
      col = colTitleX,
      font = fontTitleX),
    ylab = list(label = titleY,
      cex = cexTitleY,
      rot = rotTitleY,
      col = colTitleY,
      font = fontTitleY),
    panel = labels,
    pretty = TRUE,
    par.settings = list(panel.background = list(col = colFrame)),
    scales = list(
      x = list(cex = cexLabX,
        rot = rotLabX,
        col = colLabX,
        font = fontLabX),
      y = list(cex = cexLabY,
        rot = rotLabY,
        col = colLabY,
        font = fontLabY),
      tck = axisTicks,
      alternating = posLab),
    aspect = 'fill',
    col.regions = cols,
    cuts = 100,
    at = seq(iLowerRange, iUpperRange, 0.01),
    main = list(label = main,
      cex = cexMain,
      rot = rotMain,
      col = colMain,
      font = fontMain),
    colorkey = list(space = posColKey,
      labels = list(cex = cexLabColKey)))

  # return plot?
  if (returnPlot == TRUE) {
    return(l)
  } else if (returnPlot == FALSE) {
    l
  }
}
