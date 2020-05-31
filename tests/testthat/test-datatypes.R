skip("missing lots of things")

test_that('datatypes', {
  # pca
  expect_type(pca(mat, metadata),
    c('S4', 'list'))
  expect_type(pca(center, scale),
    c('logical'))
  expect_type(pca(removeVar),
    c('double'))
  expect_lt(pca(removeVar), 1)

  # screeplot
  expect_is(screeplot(pcaobj), class = 'pca')
  expect_type(screeplot(components, xlab, ylab, title, subtitle, caption,
    colBar, colCumulativeSumLine, colCumulativeSumPoints, hlineType,
    hlineCol, vlineType, vlineCol, borderColour),
    c('language', 'character'))
  expect_type(screeplot(drawCumulativeSumLine, drawCumulativeSumPoints,
    gridlines.major, gridlines.minor, returnPlot),
    c('logical'))
  expect_type(screeplot(xlim, ylim, xlabAngle, xlabhjust, xlabvjust,
    ylabAngle, ylabhjust, ylabvjust, axisLabSize, titleLabSize,
    subtitleLabSize, captionLabSize, sizeCumulativeSumLine,
    sizeCumulativeSumPoints, hline, hlineWidth, vline, vlineWidth,
    borderWidth),
    c('double'))

  # biplot
  expect_is(biplot(pcaobj), class = 'pca')
  expect_type(biplot(x, y, colLoadingsNames, colConnectorsLoadings,
    colLoadingsArrows, colby, colkey, singlecol, legendPosition,
    lab, selectLab, colConnectors, xlab, ylab, title, subtitle,
    caption, hlineType, hlineCol, vlineType, vlineCol, borderColour),
    c('language', 'character'))
  expect_type(biplot(showLoadings, showLoadingsNames, boxedLoadingsNames,
    drawConnectorsLoadings, boxedLabels, drawConnectors, gridlines.major,
    gridlines.minor, returnPlot),
    c('logical'))
  expect_type(biplot(ntopLoadings, sizeLoadingsNames, widthConnectorsLoadings,
    lengthLoadingsArrowsFactor, widthLoadingsArrows, alphaLoadingsArrow, shape,
    shapekey, pointSize, legendLabSize, legendIconSize, xlim, ylim, labSize,
    labhjust, labvjust, widthConnectors, xlabAngle, xlabhjust, xlabvjust,
    ylabAngle, ylabhjust, ylabvjust, axisLabSize, titleLabSize,
    subtitleLabSize, captionLabSize, hline, hlineWidth, vline, vlineWidth,
    borderWidth),
    c('double'))

  # getComponents
  expect_is(getComponents(pcaobj), class = 'pca')
  expect_type(getComponents(components),
    c('language', 'character'))

  # getVars
  expect_is(getVars(pcaobj), class = 'pca')
  expect_type(getVars(components),
    c('language', 'character'))

  # getLoadings
  expect_is(getLoadings(pcaobj), class = 'pca')
  expect_type(getLoadings(components),
    c('language', 'character'))

  # pairsplot
  expect_is(pairsplot(pcaobj), class = 'pca')
  expect_type(pairsplot(components, x, y, colby, colkey, singlecol,
    legendPosition, lab, selectLab, colConnectors, xlab, ylab,
    title, hlineType, hlineCol, vlineType, vlineCol, borderColour),
    c('language', 'character'))
  expect_type(pairsplot(triangle, plotaxeslab, drawConnectors,
    gridlines.major, gridlines.minor, returnPlot),
    c('logical'))
  expect_type(pairsplot(trianglelabSize, margingaps, ncol, nrow, shape,
    shapekey, pointSize, legendLabSize, legendIconSize, xlim, ylim,
    labSize, labhjust, labvjust, widthConnectors, xlabAngle, xlabhjust,
    xlabvjust, ylabAngle, ylabhjust, ylabvjust, axisLabSize, titleLabSize,
    hline, hlineWidth, vline, vlineWidth, borderWidth),
    c('double'))

  # plotloadings
  expect_is(plotloadings(pcaobj), class = 'pca')
  expect_type(plotloadings(components, col, legendPosition,
    positionConnectors, typeConnectors, endsConnectors, colConnectors,
    xlab, ylab, title, subtitle, caption, hlineType, hlineCol, vlineType,
    vlineCol, borderColour),
    c('language', 'character'))
  expect_type(plotloadings(absolute, drawConnectors, gridlines.major,
    gridlines.minor, returnPlot),
    c('logical'))
  expect_type(plotloadings(rangeRetain, colMidpoint, shape, shapeSizeRange,
    legendLabSize, legendIconSize, xlim, ylim, labSize, labhjust, labvjust,
    widthConnectors, lengthConnectors, xlabAngle, xlabhjust, xlabvjust,
    ylabAngle, ylabhjust, ylabvjust, axisLabSize, titleLabSize,
    subtitleLabSize, captionLabSize, hline, hlineWidth, vline, vlineWidth,
    borderWidth),
    c('double'))

  # eigencorplot
  expect_is(eigencorplot(pcaobj), class = 'pca')
  expect_type(eigencorplot(components, metavars, titleX, colTitleX, titleY,
    colTitleY, colLabX, colLabY, posLab, col, posColKey, colCorval, main,
    colMain, corFUN, corUSE, corMultipleTestCorrection, signifSymbols, colFrame),
    c('language', 'character'))
  expect_type(eigencorplot(scale, plotRsquared, returnPlot),
    c('logical'))
  expect_type(eigencorplot(cexTitleX, rotTitleX, fontTitleX, cexTitleY,
    rotTitleY, fontTitleY, cexLabX, rotLabX, fontLabX, cexLabY, rotLabY,
    fontLabY, cexLabColKey, cexCorval, fontCorval, cexMain, rotMain,
    fontMain, signifCutpoints),
    c('double'))
})
