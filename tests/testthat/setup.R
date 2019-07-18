set.seed(1000)
ngenes <- 1000
ncells <- 200
npops <- 10

mu <- matrix(rnorm(npops*ngenes), ncol=npops)
lcounts <- mu[,sample(ncol(mu), ncells, replace=TRUE)]
lcounts <- lcounts + rnorm(length(lcounts))
