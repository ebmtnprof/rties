% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rtiesShared.R
\name{histAll}
\alias{histAll}
\title{Histograms for all numeric variables in a dataframe.}
\usage{
histAll(basedata)
}
\arguments{
\item{basedata}{A user-provided dataframe.}
}
\value{
No return value. Prints plots to the console.
}
\description{
Useful for checking distributions to assess normality
}
\examples{
data <- rties_ExampleDataShort
vars <- c("reltime","ambiv","love","conflict")
newData <- data[vars ]
histAll(newData)

}
