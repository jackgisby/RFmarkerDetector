context("combinatorialCV")

data(cachexiaData)
cachexiaData <- cachexiaData[,1:10]

params <- list(ntrees = 100, nsplits = 10, test_prop = 1 / 3)
combinatorialRFMCCV_res <- combinatorialRFMCCV(expr_df = cachexiaData, 
                                               parameters = params)
combinatorialRFMCCV_res <- combinatorialRFMCCV_res[["best_model"]]

label_table <- table(as.vector(as.matrix(combinatorialRFMCCV_res$labels)))

test_that("combinatorialRFMCCV", {
    expect_equal(names(combinatorialRFMCCV_res), c("predictions", "labels", "models"))
    
    expect_equal(sum(combinatorialRFMCCV_res$predictions), 101.17)
    
    expect_equal(names(label_table), c("cachexic", "control"))
    expect_equal(as.vector(label_table), c(161, 89))
    
    expect_equal(length(combinatorialRFMCCV_res$models), 10)
    expect_equal(round(mean(combinatorialRFMCCV_res$models[[1]]$err.rate[,1]), 4),
                 0.2898)
})
