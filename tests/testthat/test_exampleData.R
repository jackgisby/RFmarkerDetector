context("exampleData")

data(simpleData)
data(cachexiaData)

test_that("data", {
    expect_equal(length(simpleData), 2)
    expect_equal(length(cachexiaData), 65)
    
    expect_equal(length(simpleData[[1]]), 200)
    expect_equal(length(cachexiaData[,1]), 77)
})
