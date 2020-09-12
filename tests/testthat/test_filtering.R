context("test_filtering")

data(cachexiaData)
cachexiaData <- cachexiaData[,1:10]

expr_df.filtered <- lqvarFilter(cachexiaData)

test_that("lqvarFilter", {
    expect_equal(length(expr_df.filtered$filtered_expr_df), 10)
    expect_equal(expr_df.filtered$filtered_expr_df[[3]][1], 40.85)
})

v <- runif(10, min = 5, max = 30)
rsd_res <- rsd(v)

test_that("rsd", {
    expect_equal(rsd_res, 34.04)
})

expr_df.filtered <- rsdFilter(cachexiaData, threshold = 15)

test_that("rsdFilter", {
    expect_equal(length(expr_df.filtered$filtered_expr_df), 10)
    expect_equal(expr_df.filtered$filtered_expr_df[[3]][1], 40.85)
})
