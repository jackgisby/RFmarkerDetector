#' Monte Carlo cross-validation of Random Forest models
#'
#' This function allows to perform a Monte Carlo cross-validation of a Random Forest
#' 
#' @param expr_df a n x p expr_dfframe used to build the models. The first two columns
#' must represent respectively the sample names and the class labels related to each sample
#' 
#' @param nsplits the number of random splittings of the original expr_df into training and test expr_df sets
#' 
#' @param test_prop the percentage (expressed as a real number) of the observations of the original expr_df
#' to be included in each test set
#' 
#' @param opt_params a list of optional parameters characterizing both the model to be validated and the input expr_df.
#' It may include parameters like the number of trees (ntree), the mtry, or the eventual reference class label (ref_level) of the expr_df
#' 
#' @return a list of three elements: \itemize{
#' \item a n x p expr_dfframe representing the predictions during the cross-validation process. Its number of lines is equal to the number of observations included in each test set
#' and the number of columns is equal to the number of test sets (defined by the nsplits input parameter).
#' \item a n x p expr_dfframe representing the labels associated to the samples of each test set. Its number of lines n is equal to the number observations in each test set while p is the number of test sets defined by the nsplits input parameter
#' \item a list of p random forest models tested in the cross-validation process
#' }
#'
#' @examples
#' data(cachexiaData)
#' params <- list(ntrees = 10, ref_level = levels(cachexiaData[, 2])[1])
#' 
#' mccv_obj <- rfMCCV(cachexiaData[,1:10], nsplits = 5, test_prop = 1 / 3, opt_params = params)
#' 
#' @author Piergiorgio Palla
#'  
#' @export

rfMCCV <- function(expr_df, nsplits, test_prop, opt_params) {
  
  ### Random Forest Default Parameters ###
  ntree <- 1000
  MTRY <- floor(sqrt(ncol(expr_df) - 2))

  if (methods::hasArg(opt_params)) {
    if (!is.null(opt_params[["ntrees"]])) {
      ntree <- opt_params$ntrees
    }

    if (!is.null(opt_params[["mtry"]])) {
      MTRY <- opt_params$mtry
    }

    if (!is.null(opt_params[["ref_level"]])) {
      ref_level <- opt_params$ref_level
      labels <- levels(expr_df[, 2])
      if (!(ref_level %in% labels)) {
        stop("A problem occurred in opt_params: check ref_level parameter")
      }
      expr_df[, 2] <- stats::relevel(expr_df[, 2], ref = ref_level)
    }
  }

  ### First of all we split the original expr_df ###
  ### Here we compute the number of observations in each test set
  ntest <- floor(test_prop * nrow(expr_df))

  ### test_index_matrix has nsplits x ntest elements. Each row contains the indexes of the observations of the original expr_df included in one of
  ### the test sets defined by the nsplits parameter.
  set.seed(1)
  test.index.matrix <- WilcoxCV::generate.split(n = nrow(expr_df), niter = nsplits, ntest = ntest)
  
  ### From this matrix we build the splits: training sets and test sets
  m <- matrix(nrow = ntest, ncol = nsplits) ### base empty matrix
  
  ### Here we initialize our output variables ####
  predictions <- data.frame(m)
  labels <- data.frame(m)
  models <- list()
  for (i in 1:nrow(test.index.matrix)) {
    indexes <- test.index.matrix[i, ]
    testset <- expr_df[indexes, ]
    trainingset <- expr_df[-indexes, ]
    model <- randomForest::randomForest(x = trainingset[, 3:ncol(trainingset)], y = trainingset[, 2], xtest = testset[, 3:ncol(testset)], ytest = testset[
      ,
      2
    ], mtry = MTRY, ntree = ntree)

    predictions[, i] <- model$test$votes[, 2]
    labels[, i] <- testset[, 2]

    models[[i]] <- model
  }

  return(mccv(predictions, labels, models))
}

#' mccv class
#'
#' A constructor function for the S3 class mccv; the mccv class encapsulates information such as predictions and abels needed to
#' plot roc curve(s) for a cross-validated random forest model
#' 
#' @param predictions a n x p expr_dfframe of the predictions collected during the cross-validation process of the random forest, with
#' n is equal to the number of samples in each test set and p is equal to the number of test sets
#' 
#' @param labels a n x p expr_dfframe of the labels of the samples included in each test set obtained splitting the original expr_df
#' during the cross-validation process of the model
#' 
#' @param models a list of p random forest models tested in the cross-validation process
#' 
#' @return an object of class mccv
#' 
#' @author Piergiorgio Palla
#' 
#' @examples
#' ## load a simple expr_df with the vectors of the predictions and the labels resulting from a CV
#' data(simpleData)
#' mccv_obj <- mccv(simpleData$predictions, simpleData$labels, models = NULL)
#' 
#' @export

mccv <- function(predictions, labels, models) {
  res <- list(predictions = predictions, labels = labels, models = models)
  class(res) <- "mccv"
  invisible(res)
}

#' Plotting single or multiple ROC curves of the cross-validated Random Forest models
#' \code{plot_mccv} allows to plot single or multiple ROC curves to characterize the performace of a cross-validated
#' Random Forest model
#' 
#' @param x an object of class mccv
#' 
#' @param y not used
#' 
#' @param opt a list containing the following optional parameters: \itemize{
#' \item avg if the mccv object represents the predictions obtained from different cross-validation runs, we can have a different roc
#' curve for each cv run. These curves can be averaged or not. Allowed values are none (plot all curves separately),
#' horizontal (horizontal averaging), vertical(vertical averaging) and threshold (threshold averaging).
#' \item colorize a logical value which indicates if the curve(s) shoud be colorized according to the cutoff.
#' }
#' 
#' @examples
#' data(cachexiaData)
#' 
#' params <- list(ntrees = 50, ref_level = levels(cachexiaData[,2])[1] )
#' mccv_obj <- rfMCCV(cachexiaData[,1:10], nsplits = 10, test_prop = 1/3, opt_params = params)
#' params = list(avg = 'vertical', colorize = FALSE)
#' 
#' mccv.plot(mccv_obj, opt = params)
#' 
#' @author Piergiorgio Palla
#' 
#' @export

mccv.plot <- function(x, y, opt = list(avg = "vertical", colorize = F)) {
  pred <- ROCR::prediction(x$predictions, x$labels)
  perf <- ROCR::performance(pred, "tpr", "fpr")

  avg_param <- "none"
  colorize <- FALSE ### Silencing R CMD CHECK NOTES

  colorize_param <- colorize
  allowed_avg <- c("none", "vertical", "horizontal", "threshold")

  if (!is.null(opt[["avg"]]) & (opt[["avg"]] %in% allowed_avg)) {
    avg_param <- opt$avg
    if (avg_param == "threshold" & !is.null(colorize)) {
      colorize_param <- opt$colorize
    }
  } else {
    stop("avg value not permitted. Allowed values: horizontal, vertical, threshold and none")
  }

  ROCR::plot(perf, avg = avg_param, colorize = colorize_param, lwd = 2, main = "ROC curves from Monte Carlo CV")
}

#' Characterizing the performance of a Random Forest model
#'
#' This function provides the accuracy and the recall of a Random Forest model
#' 
#' @param cm the confusion matrix of a validated Random Forest model
#' 
#' @return a vector contining the accuracy, the recall and the same confusion matrix of the Random Forest model
#' 
#' @author Piergiorgio Palla

forestPerformance <- function(cm) {
  # (-) (+) ctrl case --- predicted ctrl (-) cm11 | cm12 case (+) cm21 | cm22 | | actual

  ## Accuracy: Fraction of all instances correctly classified
  accuracy <- sum(diag(cm)) / sum(rowSums(cm))
  
  ## Recall (sensitivity)
  recall <- cm[2, 2] / (cm[2, 1] + cm[2, 2])

  ## The precision for a class is the number of TP (i.e. the number of items correctly labeled as belonging to that class) divided by the total
  ## number of elements labeled as belonging to that class (i.e. the sum of TP and FP, which are items incorrectly labeled as belonging to that
  ## class)
  return(c(accuracy, recall, cm))
}

#' Extracting average accuracy and recall of a list of Random Forest models
#'
#' This function provides the average accuracy and the recall of a list of Random Forest models
#' 
#' @param model_list a list of different Random Forest models
#' 
#' @return a list of the two elements: \itemize{
#' \item avg_accuracy the average accuracy
#' \item avg_recall the average recall
#' }
#' 
#' @examples
#' data(cachexiaData)
#' params <- list(ntrees = 10, ref_level = levels(cachexiaData[,2])[1])
#' 
#' mccv_obj <- rfMCCV(cachexiaData[,1:10], nsplits = 5, test_prop = 1/3, opt_params = params)
#' models <- mccv_obj$models
#' 
#' res <- rfMCCVPerf(models)
#' 
#' @author Piergiorgio Palla
#' 
#' @export

rfMCCVPerf <- function(model_list) {
  acc_vect <- c()
  recall_vect <- c()

  for (i in 1:length(model_list)) {
    if (!is.null(model_list[[i]]$test$confusion)) {
      cm <- model_list[[i]]$test$confusion
    } else {
      cm <- model_list[[i]]$confusion
    }

    tmp <- forestPerformance(cm)
    acc_vect <- append(acc_vect, tmp[1])
    recall_vect <- append(recall_vect, tmp[2])
  }

  res <- list(avg_accuracy = mean(acc_vect), avg_recall = mean(recall_vect))
}

#' Computing the average AUC
#'
#' This function computes the average AUC of the ROC curves obtained from a cross-validation procedure
#' 
#' @param predictions a vector, matrix, list, or expr_df frame containing the predictions.
#' 
#' @param labels a vector, matrix, list, or expr_df frame containing the true class labels. It must have
#' the same dimensions as predictions
#' 
#' @return the average AUC value
#' 
#' @examples
#' 
#' ## load a simple expr_df with the vectors of the predictions and the labels resulting from a CV
#' data(simpleData)
#' avg_auc <- getAvgAUC(simpleData$predictions, simpleData$labels)
#' 
#' @author Piergiorgio Palla
#' 
#' @export

getAvgAUC <- function(predictions, labels) {

  ## here we use the ROCR library to compute the AUC print(paste('pred-lab: ', predictions, labels))
  pred <- ROCR::prediction(predictions, labels)
  
  ##### AUC ####
  auc <- ROCR::performance(pred, "auc")

  ### AUC values are contained in a list enclapsulated in a slot of the object auc Here we get the slot and convert the list to an array to easily
  ### manupulate values .
  auc_values <- unlist(methods::slot(auc, "y.values"))
  
  return(round(mean(auc_values), 3))
}

#' Extracting the best performing Random Forest model
#'
#' This function allows to find the best performing Random Forest model starting from a k-combination of its input variables
#' 
#' @param combinations a k x n matrix in which n is the number of combinations of the input variables and
#' k is the size of each combination
#' 
#' @param expr_df a n x p expr_df frame of n observations and p-2 predictors. The first two columns must represent the sample
#' names and the classes associates to each sample
#' 
#' @param params a list of params useful to perform a Monte Carlo Cross validation. It should contain the following
#' expr_df: \itemize{
#' \item ntrees the number of trees of each random forest model
#' \item nsplits the number of random splittings of the original expr_df into training and test expr_df sets
#' \item test_prop the percentage (expressed as a real number) of the observations of the original expr_df
#' to be included in each test set
#' \item ref_level the assumed reference class label
#' }
#' 
#' @param ntree the number of trees of each random forest model
#' 
#' @param nsplits the number of random splittings of the original expr_df into training and test expr_df sets
#' 
#' @param test_prop the percentage (expressed as a real number) of the observations of the original expr_df
#' to be included in each test set
#' 
#' @return a list of the following elements: \itemize{
#' \item best_model_set the set of best performing Random Forest models in terms of AUC
#' \item max_auc the maximum value of AUC corresponding to those models
#' \item biomarker_set the set of metabolites (or bins) corresponding to the best performing Random Forest
#' }
#' 
#' @details The k-combinations of the input variables is represented as a k x n matrix in which k is the size of each
#' combination and n is the number of combinations of the input variables of the original expr_df.
#' Each column of the combinations matrix contains the indexes of the input variables from the original expr_df
#' The \code{getBestRFModel} extracts a expr_df from the original one considering the indexes in these columns.
#' Then it will build a Random Forest model performing a Monte Carlo CV for each expr_df.
#' The models cross-validated will be compared considering the AUC of their averaged ROC curve.
#' The function will return the best models, the maximum value of AUC and the most relevant input variables associated
#' 
#' @examples
#' data(cachexiaData)
#' 
#' expr_df <- cachexiaData[, 1:15]
#' indexes <- 3:15
#' combinations <- utils::combn(x = indexes, m = 5) # a 5 x n_of_combinations matrix
#' test_params = list(ntrees= 10, nsplits = 10, test_prop = 1/3)
#' 
#' res <- getBestRFModel(combinations, expr_df, test_params)
#' 
#' @author Piergiorgio Palla
#' 
#' @export

getBestRFModel <- function(combinations, expr_df, params, ntree=1000, nsplits=100, test_prop=1/3) {
  REF_LEVEL <- levels(expr_df[, 2])[1]

  if (methods::hasArg(params)) {
    if (!is.null(params[["ntrees"]])) {
      ntree <- params$ntrees
    }

    if (!is.null(params[["nsplits"]])) {
      nsplits <- params$nsplits
    }

    if (!is.null(params[["test_prop"]])) {
      test_prop <- params$test_prop
    }

    if (!is.null(params[["ref_level"]])) {
      REF_LEVEL <- params$ref_level
    }
  }

  expr_df[, 2] <- stats::relevel(expr_df[, 2], ref = REF_LEVEL)
  mccv_res <- list()
  auc_values <- c()
  
  ## kxn matrix with k the dimension of each combination and n the total number of k-combinations
  for (j in 1:ncol(combinations)) {

    ### selected columns
    col_idx <- combinations[, j]
    num_of_variables <- dim(combinations)[1]
    num_of_combinations <- dim(combinations)[2]

    ### extracted expr_df ##
    extracted_expr_df <- cbind(expr_df[, 1:2], expr_df[, col_idx])
    res <- rfMCCV(expr_df = extracted_expr_df, nsplits = nsplits, test_prop = test_prop, opt_params = params)
    current_auc <- getAvgAUC(res$predictions, res$labels)
    auc_values <- append(auc_values, current_auc)
    mccv_res[[j]] <- res
  }

  #### Best model params ###
  opt_metabolites <- names(expr_df)

  max_auc <- max(auc_values)
  max_idx <- which.max(auc_values)


  best_combination <- combinations[, max_idx]
  best_model_set <- mccv_res[[max_idx]]
  plotAUCvsCombinations(auc_values, num_of_variables, num_of_combinations)

  return(list(auc = max_auc, biomarker_set = best_combination, best_model_set = best_model_set))
}
