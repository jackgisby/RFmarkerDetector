context("RFMCcv")

data(cachexiaData)
cachexiaData <- cachexiaData[,1:10]

params <- list(ntrees = 10, ref_level = levels(cachexiaData[, 2])[1])
mccv_res <- rfMCCV(cachexiaData, nsplits = 5, test_prop = 1 / 3, opt_params = params)
label_table <- table(as.vector(as.matrix(mccv_res$labels)))

test_that("rfMCCV", {
    expect_equal(names(mccv_res), c("predictions", "labels", "models"))
    
    expect_equal(sum(mccv_res$predictions), 43)
    
    expect_equal(names(label_table), c("cachexic", "control"))
    expect_equal(as.vector(label_table), c(78, 47))
    
    expect_equal(length(mccv_res$models), 5)
    expect_equal(round(mean(mccv_res$models[[1]]$err.rate[,1]), 4),
                 0.3481)
})

mccv_class <- mccv(simpleData$predictions, simpleData$labels, models = NULL)
label_table <- table(as.vector(as.matrix(mccv_class$labels)))

test_that("mccv", {
    expect_equal(names(mccv_class), c("predictions", "labels", "models"))
    
    expect_equal(sum(mccv_class$predictions), 96.9344)
    
    expect_equal(names(label_table), c("0", "1"))
    expect_equal(as.vector(label_table), c(134, 66))
})

rfMCCVPerf_res <- rfMCCVPerf(mccv_res$models)

test_that("rfMCCVPerf", {
    expect_equal(round(rfMCCVPerf_res$avg_accuracy, 3), 0.661)
    expect_equal(round(rfMCCVPerf_res$avg_recall, 3), 0.497)
})

data(simpleData)
avg_auc <- getAvgAUC(simpleData$predictions, simpleData$labels)

test_that("getAvgAUC", {
    expect_equal(avg_auc, 0.464)
})

indexes <- 3:10
combinations <- combn(x = indexes, m = 5) # a 5 x n_of_combinations matrix
getBestRFModel_res <- getBestRFModel(combinations, cachexiaData, params, ntree=50, nsplits=10)
label_table <- table(as.vector(as.matrix(getBestRFModel_res$best_model_set$labels)))

test_that("getBestRFModel", {
    expect_equal(names(getBestRFModel_res$best_model_set), c("predictions", "labels", "models"))
    
    expect_equal(sum(getBestRFModel_res$best_model_set$predictions), 102.7)
    
    expect_equal(names(label_table), c("cachexic", "control"))
    expect_equal(as.vector(label_table), c(161, 89))
    
    expect_equal(length(getBestRFModel_res$best_model_set$models), 10)
    expect_equal(round(mean(getBestRFModel_res$best_model_set$models[[1]]$err.rate[,1]), 4),
                 0.4059)
    
    expect_equal(getBestRFModel_res$auc, 0.769)
    expect_equal(getBestRFModel_res$biomarker_set, c(3, 4, 8, 9, 10))
})
