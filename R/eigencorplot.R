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
  # let's make sure that anything that isnt' numeric becomes a numeric
  # find all the character columns
  # ---select all the columns which are characters
  # using base R now for dependencies sake 
    chararcter_columns = unlist(lapply(yvals, is.numeric))  
    # negate it - basically if it
    chararcter_columns = !chararcter_columns
    # select only  the names that are true 
    chararcter_columns = names(which(chararcter_columns))
    for (c in chararcter_columns){
      print(c)
      yvals[, eval(quote(c))] = as.numeric(as.factor(yvals[, eval(quote(c))]))}

    
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
