#' Principal Component Analysis
#'
#' Performs a principal component analysis based on Singular Value Decomposition, on the given data matrix
#' and returns the result as an object of the S3 class \code{pca}
#' 
#' @param X a n x p data frame of n observations and p variables.
#' 
#' @param autoscale a logical value indicating whether the variables should be autoscaled
#' 
#' @param exclude a logical value indicating whether the first two columns should be excluded from the computation.
#' The default is TRUE, because usually the first two columns of the expr_df processed represent respectively
#' the sample names and the class labels associated with the samples
#' 
#' @return an S3 object of class \code{pca} with the following components: \itemize{
#' \item scores the scores matrix
#' \item loadings the loading matrix
#' \item variances the vector of variances explained by each PC
#' \item classes the vector of the class labels associated with the samples
#' \item features the vector with the names of the input variables
#' }
#' 
#' @examples
#' data(cachexiaData)
#' pca_obj <- pca(cachexiaData[,1:10], autoscale = TRUE, exclude = TRUE)
#' 
#' @author Piergiorgio Palla
#' 
#' @export

pca <- function(X, autoscale = T, exclude = T) {
  if (exclude) {
    ### This vector contains the labels associated with the training samples ####
    X.classes <- X[, 2]

    #### Here we collect the input features
    X.features <- colnames(X)[3:dim(X)[2]]
    X <- X[, 3:dim(X)[2]]
  } else {
    X.classes <- c()
    X.features <- colnames(X)
  }

  if (autoscale) {
    X <- autoscale(X, exclude = F)
  }
  
  ### SVD factorization: X = UDV'
  X.svd <- svd(X)

  ### Loadings: V
  X.loadings <- X.svd$v
  #######

  ### Scores: (UD) - where D is a diagonal matrix
  X.scores <- X.svd$u %*% diag(X.svd$d)

  ### Here we calcultate the variance explained by each PC
  X.vars <- X.svd$d^2 / (nrow(X))
  X.totvars <- sum(X.vars)
  X.relvars <- X.vars / X.totvars
  variances <- 100 * round(X.relvars, digits = 3)

  res <- list(scores = X.scores, loadings = X.loadings, variances = variances, classes = X.classes, features = X.features)
  class(res) <- "pca"
  return(res)
}

#' PCA Scores plot
#'
#' This function creates a plot that graphically projects the original samples onto the subspce
#' spanned by the first two principal components
#' 
#' @param pca_obj an object of class pca
#' 
#' @param expr_df a n x p dataframe representing the expr_df used to create the pca object
#' 
#' @param xrange a vector of two elements indicating the range of values along the x-axis in which
#' eventually it will be specified the name of the samples
#' 
#' @param yrange a vector of two elements indicating the range of values along the y-axis in which
#' eventually it will be specified the name of the samples
#' 
#' @details The Scores plot is used for interpreting relations among the observations
#' 
#' @examples
#' data(cachexiaData)
#' pca_obj <- pca(cachexiaData[,1:10], autoscale = TRUE, exclude = TRUE)
#' pca.scores.plot(pca_obj, cachexiaData[,1:10])
#' 
#' @author Piergiorgio Palla
#' 
#' @export

pca.scores.plot <- function(pca_obj, expr_df, xrange, yrange) {
  scores <- pca_obj$scores
  variances <- pca_obj$variances
  classes <- pca_obj$classes

  plot(scores[, 1:2], type = "n", 
       xlab = paste("PC 1 (", variances[1], "%)", sep = ""), 
       ylab = paste("PC 2 (", variances[2], "%)", sep = ""))

  graphics::abline(h = 0, v = 0, col = "darkgrey")
  graphics::grid()
  
  class_labels <- levels(classes)
  classes <- as.numeric(classes)

  numOfClasses <- length(class_labels)
  color_vector <- seq(from = 10, length.out = numOfClasses)

  graphics::legend("topright", legend = class_labels, col = color_vector, pch = 17:15, cex = 0.8, pt.cex = 1.1, text.width = 2)

  # Now we create the corresponding factor
  cl <- as.factor(classes)

  cl_idx <- 1

  while (cl_idx <= length(levels(cl))) {
    class <- levels(cl)[cl_idx]
    class <- as.numeric(class)
    pt <- which(cl == class)
    if (class == 1) {
      symbol <- 17
      color <- "red"
    } else if (class == 2) {
      symbol <- 16
      color <- "green"
    } else {
      symbol <- 15
      color <- "blue"
    }
    
    graphics::points(scores[pt, 1:2], pch = symbol, col = color)
    cl_idx <- cl_idx + 1
  }

  ### Here we try to represent the meaningful points in the graph ###
  coordinates <- data.frame(scores[, 1:2])
  coordinates$samples <- expr_df$sample
  if (methods::hasArg(xrange)) {
    xcondition <- (coordinates[, 1] >= min(xrange) & coordinates[, 1] <= max(xrange))
  } else {
    xcondition <- rep(FALSE, length(coordinates[, 1]))
  }
  if (methods::hasArg(yrange)) {
    ycondition <- (coordinates[, 2] >= min(yrange) & coordinates[, 2] <= max(yrange))
  } else {
    ycondition <- rep(FALSE, length(coordinates[, 2]))
  }

  if (any(xcondition) & any(ycondition)) {
    condition <- xcondition & ycondition
  } else {
    condition <- xcondition | ycondition
  }

  if (any(condition)) {
    coordinates <- coordinates[condition, ]
    graphics::text(coordinates, labels = coordinates$sample, cex = 0.7)
  }
}

#' PCA Loadings plot
#'
#' This function plots the relation between the original variables and the subspace dimensions.
#' It is useful for interpreting relationships among variables.
#' 
#' @param pca_obj an object of class pca
#' 
#' @param nvar the number of variables to plot
#' 
#' @examples
#' data(cachexiaData)
#' pca_obj <- pca(cachexiaData[,1:10], autoscale = TRUE, exclude = TRUE)
#' pca.loadings.plot(pca_obj, nvar = 20)
#' 
#' @author Piergiorgio Palla
#' 
#' @export

pca.loadings.plot <- function(pca_obj, nvar) {
  loadings <- pca_obj$loadings
  variances <- pca_obj$variances
  classes <- pca_obj$classes
  features <- pca_obj$features

  DEFAULT_NVAR <- round(ncol(loadings) / 2)

  ### Here we set up the plot for the PCA loadings
  red_loadings <- loadings[, 1:2]
  row_modules <- sqrt(rowSums(red_loadings^2))

  ### Here we sort the vector of the modules of each row of the loading matrix and we extract the indexes
  indexes <- (sort(row_modules, decreasing = T, index.return = T))$ix
  tmp_mat <- cbind(loadings, features)

  if (methods::hasArg("nvar") && is.numeric(nvar)) {
    if (nvar <= 0 || nvar > ncol(loadings)) {
      msg <- paste("nvar is out of range: it must be greater than zero and less than or equal to", ncol(loadings))
      warning(msg)
      l_idx <- indexes[1:DEFAULT_NVAR]
    } else {
      l_idx <- indexes[1:nvar]
    }
  } else {
    l_idx <- indexes[1:DEFAULT_NVAR]
  }

  plot(loadings[l_idx, 1] * 1.5, loadings[l_idx, 2], type = "n", 
       xlab = paste("PC 1 (", variances[1], "%)", sep = ""), 
       ylab = paste("PC 2 (", variances[2], "%)", sep = ""))

  # The factor 1.2 is to create space for the text labels
  graphics::arrows(0, 0, loadings[l_idx, 1], loadings[l_idx, 2], col = "darkgray", length = 0.15, angle = 5)
  graphics::text(loadings[l_idx, 1:2], labels = features[l_idx], cex = 0.6, col = "blue")
}

#' Scree Plot
#'
#' This function creates  a graphical display of the variances against the number of principal components.
#' It is used to determine how many components should be retained in order to explain a high percentage of the variation in the data
#' 
#' @param pca_obj an object of class pca
#' 
#' @param ncomp the number of components to plot
#' 
#' @details screplot generates 2 graphs that represent respectively the relative variances
#' and the cumulative variances associated with the principal components
#' 
#' @examples
#' data(cachexiaData)
#' pca_obj <- pca(cachexiaData[,1:10], autoscale = TRUE, exclude = TRUE)
#' screeplot(pca_obj, ncomp = 10)
#' 
#' @author Piergiorgio Palla
#' 
#' @export

screeplot <- function(pca_obj, ncomp) {
  var <- pca_obj$variances
  DEFAULT_NCOMP <- round(length(var) / 2)


  if (methods::hasArg("ncomp") && is.numeric(ncomp)) {
    if (ncomp <= 0 || ncomp > length(var)) {
      print("OK")
      msg <- paste("ncomp is out of range: it must be within range (0:", length(var), ")", sep = "")
      warning(msg)
      ncomp <- DEFAULT_NCOMP
    } else {
      ncomp <- ncomp
    }
  } else {
    ncomp <- DEFAULT_NCOMP
  }
  
  graphics::par(mfrow = c(2, 1))
  max_var <- var[1]

  var <- var[1:ncomp]
  bplt1 <- graphics::barplot(var, main = "Relative variances", names.arg = paste("PC", 1:ncomp), ylim = c(0, 1.5 * max_var), col = 1:10)

  graphics::text(bplt1, y = var + 2.5, labels = paste(var, "%"), xpd = T, cex = 0.6)

  bplt2 <- graphics::barplot(cumsum(var), main = "Cumulative variances (%)", names.arg = paste("PC", 1:ncomp), ylim = c(0, 100), col = 10:1)
  graphics::text(bplt2, y = cumsum(var) + 5, labels = paste(cumsum(var), "%"), xpd = T, cex = 0.6)
}

#' mds class
#'
#' A constructor function for the S3 class mds; the mds class encapsulates useful information for generating MDS plots
#' 
#' @param expr_df a n x p dataframe representing the expr_df to be used to train the model. The first two columns
#' must represent respectively the samples names and the class labels related to each sample
#' 
#' @param opt a list of optional parameters useful to train the random forest. It may include the number
#' of trees (ntree), the parameter mtry and the seed
#' 
#' @return an object of class mds including the following attributes: \itemize{
#' \item model an object of class random forest containing the proximity matrix
#' \item classes a factor with the classes associated with each sample, used to train the model
#' \item sample_names the vector of the names of the samples
#' }
#' 
#' @author Piergiorgio Palla
#' 
#' @examples
#' data(cachexiaData)
#' params <- list(ntree = 50, mtry = round(sqrt(ncol(cachexiaData[,1:10]) - 2)), seed = 1)
#' mds_obj <- mds(cachexiaData[,1:10], opt = params)
#' 
#' @export

mds <- function(expr_df, opt = list(ntree = 1000, mtry = round(sqrt(ncol(expr_df) - 2)), seed = 1)) {
  set.seed(opt$seed)
  model <- randomForest::randomForest(expr_df[, 3:ncol(expr_df)], expr_df[, 2], mtry = opt$mtry, ntree = opt$ntree, proximity = T)
  classes <- expr_df[, 2]
  sample_names <- expr_df[, 1]

  res <- list(model = model, labels = classes, sample_names = sample_names)
  class(res) <- "mds"
  invisible(res)
}

#' Multi-dimensional Scaling (MDS) Plot
#'
#' This function plots the scaling coordinates of the proximity matrix from random forest.
#' 
#' @param mds_obj an object of class mds
#' 
#' @param xrange a vector of two elements indicating the interval along the x-axis in which we want to display the names of the samples
#' 
#' @param yrange a vector of two elements indicating the interval along the y-axis in which we want to display the names of the samples
#' 
#' @details From the trained model we can get the dissimilarity matrix  ** 1 - prox(i,j) **
#' The entries of this matrix can be seen as squared distances in a Euclidean high dimensional space.
#' After having calculated scaling coordinates, we can project the data onto a lower dimensional space, preserving
#' (as much as possible) the distances between the orginal points.
#' This plot can be useful for discovering patterns in data.
#' 
#' @examples
#' data(cachexiaData)
#' params = list(ntree = 50, mtry = round(sqrt(ncol(cachexiaData) - 2)), seed = 1)
#' mds_obj <- mds(cachexiaData[,1:10], opt = params)
#' mds.plot(mds_obj = mds_obj)
#' 
#' @export

mds.plot <- function(mds_obj, xrange, yrange) {
  graphics::par(mfrow = c(1, 1))
  
  ### DMat is the dissimilarity matrix
  model <- mds_obj$model
  samples <- mds_obj$sample_names
  classes <- mds_obj$labels
  print(classes)
  
  n_of_classes <- length(levels(classes))

  colors <- c("red", "blue", "darkgreen", "brown", "darkmagenta", "deeppink", "maroon", "orange", "turquoise", "violet")

  DMat <- randomForest::MDSplot(model, classes,
    palette = colors[1:n_of_classes], pch = as.numeric(classes), xlab = "1st scaling coordinate",
    ylab = "2nd scaling coordinate"
  )

  class_labels <- levels(classes)

  graphics::legend("topleft", legend = class_labels, col = colors[1:n_of_classes], pch = seq_along(levels(classes)), cex = 0.4)
  graphics::title(main = "Multi Dimensional Scaling Plot of RF Proximity Matrix")

  coordinates <- data.frame(DMat$points)
  coordinates$samples <- samples

  # Here we test if the arguments xrange and yrange has been provided
  if (methods::hasArg(xrange)) {
    xcondition <- (coordinates[, 1] >= min(xrange) & coordinates[, 1] <= max(xrange))
  } else {
    xcondition <- rep(FALSE, length(coordinates[, 1]))
  }
  if (methods::hasArg(yrange)) {
    ycondition <- (coordinates[, 2] >= min(yrange) & coordinates[, 2] <= max(yrange))
  } else {
    ycondition <- rep(FALSE, length(coordinates[, 2]))
  }

  if (any(xcondition) & any(ycondition)) {
    condition <- xcondition & ycondition
  } else {
    condition <- xcondition | ycondition
  }

  if (any(condition)) {
    coordinates <- coordinates[condition, ]
    graphics::text(coordinates, labels = coordinates$sample, cex = 0.7)
  }
}
