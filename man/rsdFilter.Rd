% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/filtering.R
\name{rsdFilter}
\alias{rsdFilter}
\title{Filtering less informative variables}
\usage{
rsdFilter(expr_df, threshold, exclude = T)
}
\arguments{
\item{expr_df}{a n x p expr_df frame with n observations and p columns. While the first two columns usually represent the names of the samples and the
class labels related to each sample respectively, the remaining columns represent metabolite concentrations measured by 1H NMR or bins of 1H NMR spectra}

\item{threshold}{a numeric value representing a limit: each predictor with a relative standard deviation lower than that
will be removed form the dataframe}

\item{exclude}{a logical variable which stores a simple True / False setting. If set to True the filtering method will exclude the first two columns.}
}
\value{
a list containing the filtered expr_df, a vector with the names of the varables excluded and a vector with the indexes of the variables eliminated
}
\description{
\code{rsdFilter} removes from the expr_dfframe the predictors with a relative standard deviation less
than or equal to an inserted threshold
}
\examples{

## load the included example data
data(cachexiaData)

## call rsdFilter with the parameter exclude set to TRUE (default)
## in order to exclude the first two columns of the expr_df from scaling
expr_df.filtered <- rsdFilter(cachexiaData[,1:10], threshold = 15, exclude = TRUE)

}
\author{
Piergiorgio Palla
}
