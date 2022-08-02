

######### The "biVar" file includes functions that are focused on modeling dynamics between two different variables within-person (e.g., bivariate data), rather than the same variable between two partners in a dyad (e.g., dyadic data, which is the focus of most of rties)

############## makeBiVarData

#' Takes typical long-format time-nested-in-person data, and stacks two user-chosen observed variables on top of each other so they can be treated as "bivariate" within person. In other words, two time-series variables from each person are stacked on top of each other, forming a bivariate pair of variables within person (e.g., time in variable in person). 
#' 
#' The resulting data can either: 1) be used with many of the other rties functions (a number of the preparatory functions and plotting will work in this way), but instead of "dyad" being the highest nesting variable, "person" is and should be substituted instead of dyad wherever you would otherwise use dyad, or 2) the data can be reformated again with the "dataPrep" function in rties (see overview_data_prep vignette for information on how to do this), with the resulting data in a format that can be used with any of the other rties functions. 
#'
#' @param basedata The original dataframe provided by the user that includes at least two time-series variables nested within-person
#' @param personId The name of the column in the dataframe that has the person-level identifier.
#' @param time_name The name of the column in the dataframe that indicates sequential temporal observations.
#' @param obs1_name The name of the column in the dataframe that has the first time-series variable to be stacked.
#' #' @param obs2_name The name of the column in the dataframe that has the second time-series variable to be stacked.
#' @param labels A string vector with the names of the variables that are being stacked.
#' @param idConvention A value that will be added to the varId of the first variable to get the  varId for the second variable. It should be a larger value than the highest personId. This varId will then be used by rties in a way similar to personId when partners are nested in dyads.
#'
#' @return A dataframe that contains the original data, plus the following columns: 1) var: the names of the stacked variables (taken from "labels"). 2) obs: the stacked observed variable scores, 3) dist: a zero/one distinguishing variable, and 4) varId: a variable ID that is similar to personId for use with rties. The varId for the first stacked variable is the same as the personId, with the varId for the second stacked variable being personId + idConvention.

#' @export

makeBiVarData <- function(basedata, personId, time_name, obs1_name, obs2_name, labels, idConvention)
{
  colnames(basedata)[colnames(basedata)== personId] <- "person"
  colnames(basedata)[colnames(basedata)== time_name] <- "time"
  colnames(basedata)[colnames(basedata)== obs1_name] <- "obs1"
  colnames(basedata)[colnames(basedata)== obs2_name] <- "obs2"
  
  newData <- reshape(basedata, varying=c("obs1","obs2"), timevar="var", idvar=c("person","time"), direction="long", sep="")
  
  ### the next 3 steps must be in order
  # make a zero/one distinguishing variable
  newData$dist <- newData$var - 1
  
  # make a variable ID, nested within person, that will function like a person ID in rties code
  newData$varId <- ifelse(newData$dist == 0, newData$person, newData$person + idConvention)
  
  # make var into a factor with recognizable names
  newData$var <- factor(newData$var, levels=c(1,2), labels=labels)
  
  colnames(newData)[colnames(newData)== "person"] <- personId
  colnames(newData)[colnames(newData)== "time"] <- time_name
  colnames(newData)[colnames(newData)== "obs1"] <- obs1_name
  colnames(newData)[colnames(newData)== "obs12"] <- obs2_name
  
  return(newData)
}



############## makeCrossCorBiVar

#' Takes typical time-series wide-format data (e.g., multiple time-varying variables for each person in wide format) and calculates cross-correlations for two user-specified variables within a specified maximum number of lags. It returns a dataframe with the largest absolute cross-correlation and its lag added for each person (e.g., it returns either the most negative or most positive cross-correlation, whichever is larger in absolute terms -- the sign is retained).
#'
#' @param basedata The original dataframe provided by the user that includes at least two time-series variables nested within-person
#' @param personId The name of the column in the dataframe that has the person-level identifier.
#' @param time_name The name of the column in the dataframe that indicates sequential temporal observations.
#' @param obs1_name The name of the column in the dataframe that has the first time-series variable to be stacked.
#' @param obs2_name The name of the column in the dataframe that has the second time-series variable to be stacked.
#'
#' @return A cross-sectional version of the original dataframe with maximal absolute-value cross-correlations and their lags added for each person.

#' @export

makeCrossCorBiVar <- function(basedata, personId, obs1_name, obs2_name, maxLag){
  
  colnames(basedata)[colnames(basedata)== personId] <- "person"
  colnames(basedata)[colnames(basedata)== obs1_name] <- "a"
  colnames(basedata)[colnames(basedata)== obs2_name] <- "b"
  
  newId <- unique(factor(basedata$person))
  maxCC <- list()
  
  for(i in 1:length(newId)){
    datai <- basedata[basedata$person == newId[i],]
    a <- datai$a
    b <- datai$b
    maxCC[[i]] <- signAbsMaxCC(a, b, maxLag = maxLag)
  }		
  
  maxCC <- as.data.frame(do.call(rbind, maxCC)) 	
}

################ signAbsMaxCC

#' A helper function for makeCrossCorBiVar

#' @param a First time-series used in the cross-correlation
#' @param b Second time-series used in the cross-correlation
#' @param lagMax Maximum lag at which to calculate the acf. Default is 10*log10(N/m) where N is the number of observations and m the number of series. 
#' 
#' @return A list of maximum absolute value cross-correlations and the lag at which they occurred.

#' @importFrom stats ccf na.exclude

signAbsMaxCC <- function (a, b, maxLag) {
  d <- ccf(a, b, plot = FALSE, na.action=na.exclude, lag.max = maxLag)
  cor = d$acf[ ,,1] 
  lag = d$lag[ ,,1] 
  corLag = data.frame(cor,lag) 
  for(j in 1:nrow(corLag)) if(is.na(corLag[j,1])) results=NA else {
    resultMax <- corLag[which.max(corLag$cor),] 
    resultMin <- corLag[which.min(corLag$cor),]
    resultMinAbs <- abs(resultMin)
    resultMaxAbs <- abs(resultMax)
    resultAbs = which.max(c(resultMinAbs$cor, resultMaxAbs$cor))
    if (resultAbs ==1) results =resultMin else results = resultMax}
  output <-data.frame(results = results)
  output
}










