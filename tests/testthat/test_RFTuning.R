context("test_RFTuning")

data(cachexiaData)
cachexiaData <- cachexiaData[,1:10]

tuneMTRY_res <- tuneMTRY(cachexiaData, iterations = 10, maxntree = 50, mtry_length = 5, graph = FALSE)

test_that("tuneMTRY", {
    expect_equal(round(mean(tuneMTRY_res$oob), 3), 0.327)
})

optimizeMTRY_res <- optimizeMTRY(tuneMTRY_res$oob)

test_that("optimizeMTRY", {
  expect_equal(length(optimizeMTRY_res$mean_matrix), 5)
  expect_equal(length(optimizeMTRY_res$ci_matrix[1,]), 5)
  expect_equal(length(optimizeMTRY_res$sd_matrix), 5)
})

tuneNTREE_res <- tuneNTREE(cachexiaData, 3, 5, minNTREE = 10, pace = 10, seq_length = 5)

test_that("tuneNTREE", {
  expect_equal(round(mean(tuneNTREE_res), 3), 0.326)
})
