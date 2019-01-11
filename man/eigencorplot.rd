\name{eigencorplot}

\alias{eigencorplot}

\title{eigencorplot}

\description{Correlate principal components to continuous variable metadata and test significancies of these.}

\usage{
eigencorplot(pcaobj,
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
  signifSymbols = c('***', '**', '*', ''),
  signifCutpoints = c(0, 0.001, 0.01, 0.05, 1),
  colFrame = 'white',
  plotRsquared = FALSE,
  returnPlot = TRUE)
}

\arguments{
  \item{pcaobj}{Object of class 'pca' created by pca(). REQUIRED.}
  \item{components}{The principal components to be included in the plot.
  DEFAULT = getComponents(pcaobj, seq_len(10)). OPTIONAL.}
  \item{metavars}{A vector of column names in metadata representing continuos
  variables. REQUIRED.}
  \item{titleX}{X-axis title. DEFAULT = ''. OPTIONAL.}
  \item{cexTitleX}{X-axis title cex. DEFAULT = 1.0. OPTIONAL.}
  \item{rotTitleX}{X-axis title rotation in degrees. DEFAULT = 0. OPTIONAL.}
  \item{colTitleX}{X-axis title colour. DEFAULT = 'black'. OPTIONAL.}
  \item{fontTitleX}{X-axis title font style. 1, plain; 2, bold; 3, italic; 4,
  bold-italic. DEFAULT = 2. OPTIONAL.}
  \item{titleY}{Y-axis title. DEFAULT = ''. OPTIONAL.}
  \item{cexTitleY}{Y-axis title cex. DEFAULT = 1.0. OPTIONAL.}
  \item{rotTitleY}{Y-axis title rotation in degrees. DEFAULT = 0. OPTIONAL.}
  \item{colTitleY}{Y-axis title colour. DEFAULT = 'black'. OPTIONAL.}
  \item{fontTitleY}{Y-axis title font style. 1, plain; 2, bold; 3, italic; 4,
  bold-italic. DEFAULT = 2. OPTIONAL.}
  \item{cexLabX}{X-axis labels cex. DEFAULT = 1.0. OPTIONAL.}
  \item{rotLabX}{X-axis labels rotation in degrees. DEFAULT = 0. OPTIONAL.}
  \item{colLabX}{X-axis labels colour. DEFAULT = 'black'. OPTIONAL.}
  \item{fontLabX}{X-axis labels font style. 1, plain; 2, bold; 3, italic; 4,
  bold-italic. DEFAULT = 2. OPTIONAL.}
  \item{cexLabY}{Y-axis labels cex. DEFAULT = 1.0. OPTIONAL.}
  \item{rotLabY}{Y-axis labels rotation in degrees. DEFAULT = 0. OPTIONAL.}
  \item{colLabY}{Y-axis labels colour. DEFAULT = 'black'. OPTIONAL.}
  \item{fontLabY}{Y-axis labels font style. 1, plain; 2, bold; 3, italic; 4,
  bold-italic. DEFAULT = 2. OPTIONAL.}
  \item{posLab}{Positioning of the X- and Y-axis labels. 'bottomleft', bottom
  and left; 'topright', top and right; 'all', bottom / top and left /right;
  'none', no labels. DEFAULT = 'bottomleft'. OPTIONAL.}
  \item{col}{Colour shade gradient for RColorBrewer. DEFAULT = c('blue4',
  'blue3', 'blue2', 'blue1', 'white', 'red1', 'red2', 'red3', 'red4').
  OPTIONAL.}
  \item{posColKey}{Position of colour key. 'bottom', 'left', 'top', 'right'.
  DEFAULT = 'right'. OPTIONAL.}
  \item{cexLabColKey}{Colour key labels cex. DEFAULT = 1.0. OPTIONAL.}
  \item{cexCorval}{Correlation values cex. DEFAULT = 1.0. OPTIONAL.}
  \item{colCorval}{Correlation values colour. DEFAULT = 'black'. OPTIONAL.}
  \item{fontCorval}{Correlation values font style. 1, plain; 2, bold; 3,
  italic; 4, bold-italic. DEFAULT = 1. OPTIONAL.}
  \item{scale}{Logical, indicating whether or not to scale the colour range
  to max and min cor values. DEFAULT = TRUE. OPTIONAL.}
  \item{main}{Plot title. DEFAULT = ''. OPTIONAL.}
  \item{cexMain}{Plot title cex. DEFAULT = 2. OPTIONAL.}
  \item{rotMain}{Plot title rotation in degrees. DEFAULT = 0. OPTIONAL.}
  \item{colMain}{Plot title colour. DEFAULT = 'black'. OPTIONAL.}
  \item{fontMain}{Plot title font style. 1, plain; 2, bold; 3, italic; 4,
  bold-italic. DEFAULT = 2. OPTIONAL.}
  \item{corFUN}{Correlation method: 'pearson', 'spearman', or 'kendall'.
  DEFAULT = 'pearson'. OPTIONAL.}
  \item{corUSE}{Method for handling missing values (see documentation for cor
  function via ?cor). 'everything', 'all.obs', 'complete.obs',
  'na.or.complete', or 'pairwise.complete.obs'. DEFAULT =
  'pairwise.complete.obs'. OPTIONAL.}
  \item{signifSymbols}{Statistical significance symbols to display beside
  correlation values. DEFAULT = c('***', '**', '*', ''). OPTIONAL.}
  \item{signifCutpoints}{Cut-points for statistical significance. DEFAULT =
  c(0, 0.001, 0.01, 0.05, 1). OPTIONAL.}
  \item{colFrame}{Frame colour. DEFAULT = 'white'. OPTIONAL.}
  \item{plotRsquared}{Logical, indicating whether or not to plot R-squared
  values. DEFAULT = FALSE. OPTIONAL.}
  \item{returnPlot}{Logical, indicating whether or not to return the plot
  object. DEFAULT = TRUE. OPTIONAL.}
}

\value{
A \code{\link{lattice}} object.
}

\author{
Kevin Blighe <kevin@clinicalbioinformatics.co.uk>, Myles Lewis <myles.lewis@qmul.ac.uk>
}

\examples{
  options(scipen=10)
  options(digits=6)

  col <- 20
  row <- 20000
  mat1 <- matrix(
    rexp(col*row, rate = 0.1),
    ncol = col)
  rownames(mat1) <- paste0('gene', 1:nrow(mat1))
  colnames(mat1) <- paste0('sample', 1:ncol(mat1))

  mat2 <- matrix(
    rexp(col*row, rate = 0.1),
    ncol = col)
  rownames(mat2) <- paste0('gene', 1:nrow(mat2))
  colnames(mat2) <- paste0('sample', (ncol(mat1)+1):(ncol(mat1)+ncol(mat2)))

  mat <- cbind(mat1, mat2)

  metadata <- data.frame(row.names = colnames(mat))
  metadata$Group <- rep(NA, ncol(mat))
  metadata$Group[seq(1,40,2)] <- 'A'
  metadata$Group[seq(2,40,2)] <- 'B'
  metadata$CRP <- sample.int(100, size=ncol(mat), replace=TRUE)
  metadata$ESR <- sample.int(100, size=ncol(mat), replace=TRUE)

  p <- pca(mat, metadata = metadata, removeVar = 0.1)

  eigencorplot(p, components = getComponents(p, 1:10),
    metavars = c('ESR', 'CRP'))
}
