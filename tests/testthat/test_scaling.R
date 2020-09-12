context("scaling")

data(cachexiaData)
cachexiaData <- cachexiaData[,1:10]

par_res <- paretoscale(cachexiaData)
auto_res <- autoscale(cachexiaData)
mean_res <- meanCenter(cachexiaData)

test_that("scaling", {
    expect_equal(nrow(par_res), 77)
    expect_equal(ncol(par_res), 10)
    expect_equal(round(par_res$X1.6.Anhydro.beta.D.glucose[1], 2), -5.68)
    
    expect_equal(nrow(auto_res), 77)
    expect_equal(ncol(auto_res), 10)
    expect_equal(round(auto_res$X1.6.Anhydro.beta.D.glucose[1], 2), -0.5)
    
    expect_equal(nrow(mean_res), 77)
    expect_equal(ncol(mean_res), 10)
    expect_equal(round(mean_res$X1.6.Anhydro.beta.D.glucose[1], 2), -64.78)
})
