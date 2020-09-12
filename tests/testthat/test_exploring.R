context("exploring")

data(cachexiaData)
cachexiaData <- cachexiaData[,1:10]

pca_res <- pca(cachexiaData)

params <- list(ntree = 10, mtry = round(sqrt(ncol(cachexiaData) - 2)), seed = 1)
mds_obj <- mds(cachexiaData, opt = params)

test_that("exploratory_plots", {
    expect_equal(length(pca_res), 5)
    expect_equal(pca_res$scores[[1]], 1.92064)
    
    expect_equal(mds_obj$model$err.rate[[1]], 0.36)
})
