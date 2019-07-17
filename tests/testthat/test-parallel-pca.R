# This checks the parallelPCA function.
# require(PCAtools); require(testthat); source("test-parallel-pca.R")

set.seed(1001)
test_that("parallelPCA works as expected", {
    threshold <- 0.1
    pcs <- parallelPCA(lcounts, niters=20, threshold=threshold)
    permuted <- pcs$permuted
    original <- pcs$original$variance

    pvals <- rowMeans(permuted >= original)
    expect_identical(pcs$n, min(which(pvals > threshold)-1L))
    expect_true(sum(original) <= 100)
    expect_true(all(rowSums(permuted) <= 100))
    expect_identical(ncol(permuted), 20L)

    var.exp <- prcomp(t(lcounts))$sdev^2
    frac.exp <- var.exp/sum(var.exp)
    expect_equal(frac.exp[seq_along(original)], unname(original)/100)

    # Respects max.rank when choosing the number of PCs.
    set.seed(100)
    pcs.x <- parallelPCA(lcounts, niters=3, max.rank=5)
    expect_identical(ncol(pcs.x$original$rotated), 5L)
    expect_identical(pcs.x$original$rotated, pcs$original$rotated[,1:5])
})

test_that("parallelPCA respects the seed", {
    set.seed(100)
    pcs <- parallelPCA(lcounts, niters=3)
    set.seed(100)
    pcs2 <- parallelPCA(lcounts, niters=3)
    expect_identical(pcs2, pcs)
    pcs3 <- parallelPCA(lcounts, niters=3)
    expect_false(identical(pcs3, pcs))

    # With irlba:
    set.seed(100)
    ipcs <- parallelPCA(lcounts, niters=3, BSPARAM=BiocSingular::IrlbaParam(), max.rank=10)
    set.seed(100)
    ipcs2 <- parallelPCA(lcounts, niters=3, BSPARAM=BiocSingular::IrlbaParam(), max.rank=10)
    expect_identical(pcs, pcs2)
    expect_identical(ncol(ipcs), ncol(pcs))

    # With parallelization.
    BPPARAM <- BiocParallel::SnowParam(3) # define BEFORE set.seed, otherwise this sets its own seed.
    set.seed(100)
    alt <- parallelPCA(lcounts, niters=3, BPPARAM=BPPARAM)
    expect_identical(alt, pcs)
})

set.seed(1002)
test_that("parallelPCA's C++ code works as expected", {
    trans <- t(lcounts)
    shuffled <- PCAtools:::shuffle_matrix(trans, 1, 1L)
    expect_false(identical(shuffled, trans))
    expect_identical(apply(shuffled, 2, sort), apply(trans, 2, sort))

    # Is reproducible.
    shuffled2 <- PCAtools:::shuffle_matrix(trans, 1, 1L)
    expect_identical(shuffled, shuffled2)

    # Responds to the seed.
    shuffled3 <- PCAtools:::shuffle_matrix(trans, 2, 1L)
    expect_false(identical(shuffled, shuffled3))
    expect_identical(apply(shuffled, 2, sort), apply(shuffled3, 2, sort))

    # Responds to the stream.
    shuffled4 <- PCAtools:::shuffle_matrix(trans, 1, 2L)
    expect_false(identical(shuffled, shuffled4))
    expect_identical(apply(shuffled, 2, sort), apply(shuffled4, 2, sort))
})

set.seed(1003)
test_that("parallelPCA responds to different settings", {
    set.seed(100)
    pcs <- parallelPCA(lcounts, niters=3)
    set.seed(100)
    pcs1 <- parallelPCA(t(lcounts), transposed=TRUE, niters=3)
    expect_identical(pcs, pcs1)

    set.seed(100)
    pcs2 <- parallelPCA(lcounts, scale=TRUE, niters=3)
    expect_false(identical(pcs$original, pcs2$original))
    expect_false(identical(pcs$permuted, pcs2$permuted))

    set.seed(100)
    pcs3 <- parallelPCA(lcounts, removeVar=0.5, niters=3)
    expect_false(identical(pcs$original, pcs3$original))
    expect_false(identical(pcs$permuted, pcs3$permuted))
})
