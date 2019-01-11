\name{pca}

\alias{pca}

\title{pca}

\description{Principal Components Analysis (PCA) is a very powerful technique that has wide applicability in data science, bioinformatics, and further afield. It was initially developed to analyse large volumes of data in order to tease out the differences/relationships between the logical entities being analysed. It extracts the fundamental structure of the data without the need to build any model to represent it. This 'summary' of the data is arrived at through a process of reduction that can transform the large number of variables into a lesser number that are uncorrelated, i.e., the principal components', whilst at the same time being capable of easy interpretation on the original data.}

\usage{
pca(
  mat,
  metadata = NULL,
  center = TRUE,
  scale = FALSE,
  removeVar = NULL)
}

\arguments{
  \item{mat}{A data-matrix or data-frame containing numerical data only.
  REQUIRED.}
  \item{metadata}{A data-matrix or data-frame containing metadata. This will
  be stored in the resulting pca object. Strictly enforced that
  rownames(metadata) == colnames(mat). DEFAULT = NULL. OPTIONAL.}
  \item{center}{Center the data before performing PCA? Same as prcomp()
  'center' parameter. DEFAULT = TRUE. OPTIONAL.}
  \item{scale}{Scale the data? Same as prcomp() 'scale' parameter. DEFAULT
  = FALSE. OPTIONAL.}
  \item{removeVar}{Remove this % of variables based on low variance.
  DEFAULT = NULL. OPTIONAL.}
}

\value{
A \code{\link{pca}} object.
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

  getComponents(p)

  getVars(p)

  getLoadings(p)

  screeplot(p)

  screeplot(p, hline = 80)

  biplot(p)

  biplot(p, colby = 'Group', shape = 'Group')

  biplot(p, colby = 'Group', colkey = c(A = 'forestgreen', B = 'gold'),
    legendPosition = 'right')

  biplot(p, colby = 'Group', colkey = c(A='forestgreen', B='gold'),
    shape = 'Group', shapekey = c(A=10, B=21), legendPosition = 'bottom')

  pairsplot(p, triangle = TRUE)

  plotloadings(p, drawConnectors=TRUE)

  eigencorplot(p, components = getComponents(p, 1:10),
    metavars = c('ESR', 'CRP'))
}
