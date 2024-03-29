% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rtiesShared.R
\name{makeFullData}
\alias{makeFullData}
\title{Combines profile membership data from the latent profile analysis with other data for using the profile membership to predict and be predicted by the system variable.}
\usage{
makeFullData(basedata, dyadId, personId, dist_name, lpaData)
}
\arguments{
\item{basedata}{The original dataframe provided by the user that includes all variables needed for an rties analysis, including potential system and control variables, etc.}

\item{dyadId}{The name of the column in the dataframe that has the couple-level identifier.}

\item{personId}{The name of the column in the dataframe that has the person-level identifier.}

\item{dist_name}{The name of the column in the dataframe that has a variable that distinguishes the partners (e.g., sex, mother/daughter, etc) that is numeric and scored 0/1.}

\item{lpaData}{The object returned by the "inspectProfiles" function}
}
\value{
A dataframe that contains all variables needed for using the profiles to predict, or be predicted by, the system variable.
}
\description{
Combines profile membership data from the latent profile analysis with other data for using the profile membership to predict and be predicted by the system variable.
}
\examples{
data <- rties_ExampleDataShort
newData <- dataPrep(basedata=data, dyadId="couple", personId="person", 
obs_name="dial", dist_name="female", time_name="time", time_lag=2)
ic <- indivInertCoord(prepData=newData, whichModel="inertCoord")
profiles <- inspectProfiles(whichModel="inertCoord", prepData=newData, 
paramEst=ic$params, n_profiles=2)
fullData <- makeFullData(basedata=data, dyadId="couple", personId="person", 
dist_name="female", lpaData=profiles)
head(fullData)

}
