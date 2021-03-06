% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/filtering.R
\name{lqvarFilter}
\alias{lqvarFilter}
\title{Filtering 'low quality' variables from the original expr_df}
\usage{
lqvarFilter(expr_df, threshold = 0.5, exclude = T)
}
\arguments{
\item{expr_df}{a n x p expr_df frame with n observations and p columns. While the first two columns usually represent the names of the samples and the
class labels related to each sample respectively, the remaining columns represent metabolite concentrations measured by 1H NMR or bins of 1H NMR spectra}

\item{threshold}{the percentage of zero values of a variable above which it will be eliminated from the expr_df (default: 0.50)}

\item{exclude}{a logical variable which stores a simple True / False setting. If set to True the filtering method will exclude the first two columns.}
}
\value{
a list containing the filtered expr_df, a vector with the names of the varables excluded and a vector with the indexes of the variables eliminated
}
\description{
This function takes the original expr_df and filters those variables with a definite (usually relevant) percentage of zero-values
}
\examples{

## load the included example expr_df
data(cachexiaData)

## call lqvarFilter with the parameter exclude set to TRUE (default)
## in order to exclude the first two columns of the expr_df from scaling
res <- lqvarFilter(cachexiaData[,1:10], threshold = 0.4, exclude = TRUE)
expr_df.filtered <- res$filtered_expr_df

}
\author{
Piergiorgio Palla
}
