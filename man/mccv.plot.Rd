% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RFMCcv.R
\name{mccv.plot}
\alias{mccv.plot}
\title{Plotting single or multiple ROC curves of the cross-validated Random Forest models
\code{plot_mccv} allows to plot single or multiple ROC curves to characterize the performace of a cross-validated
Random Forest model}
\usage{
mccv.plot(x, y, opt = list(avg = "vertical", colorize = F))
}
\arguments{
\item{x}{an object of class mccv}

\item{y}{not used}

\item{opt}{a list containing the following optional parameters: \itemize{
\item avg if the mccv object represents the predictions obtained from different cross-validation runs, we can have a different roc
curve for each cv run. These curves can be averaged or not. Allowed values are none (plot all curves separately),
horizontal (horizontal averaging), vertical(vertical averaging) and threshold (threshold averaging).
\item colorize a logical value which indicates if the curve(s) shoud be colorized according to the cutoff.
}}
}
\description{
Plotting single or multiple ROC curves of the cross-validated Random Forest models
\code{plot_mccv} allows to plot single or multiple ROC curves to characterize the performace of a cross-validated
Random Forest model
}
\examples{
data(cachexiaData)

params <- list(ntrees = 50, ref_level = levels(cachexiaData[,2])[1] )
mccv_obj <- rfMCCV(cachexiaData[,1:10], nsplits = 10, test_prop = 1/3, opt_params = params)
params = list(avg = 'vertical', colorize = FALSE)

mccv.plot(mccv_obj, opt = params)

}
\author{
Piergiorgio Palla
}
