# Tests the pca() function.
# library(testthat); library(PCAtools); source("setup.R"); source("test-pca.R")

test_that("pca settings work as expected", {
    basic <- pca(lcounts)
    expect_equal(sum(basic$variance), 100)
    expect_equal(as.matrix(basic$rotated), prcomp(t(lcounts))$x)

    scaled <- pca(lcounts, scale=TRUE)
    expect_equal(as.matrix(scaled$rotated), prcomp(t(lcounts), scale=TRUE)$x)
    
    trunc <- pca(lcounts, rank=10)
    expect_equal(as.matrix(trunc$rotated), prcomp(t(lcounts), rank=10)$x)

    expect_identical(pca(t(lcounts), transposed=TRUE), basic)
})

test_that("removal of low-variance genes works", {
    keep <- order(DelayedMatrixStats::rowVars(lcounts), decreasing=TRUE)[seq_len(nrow(lcounts)/2)]
    ref <- pca(lcounts[keep,])
    alt <- pca(lcounts, removeVar=0.5)
    expect_equal(ref$variance, alt$variance)
})

test_that("pca works with alternative SVD algorithms", {
    basic <- pca(lcounts, rank=5)

    set.seed(0)
    irlba <- pca(lcounts, rank=5, BSPARAM=BiocSingular::IrlbaParam())
    expect_false(identical(irlba, basic))
    expect_equal(irlba$variance, basic$variance)
})

test_that("pca works with alternative matrix representations", {
    library(Matrix)
    X <- as(lcounts, "dgeMatrix")
    ref <- pca(lcounts, rank=5)
    basic <- pca(X, rank=5)
    expect_equal(ref, basic)
})
