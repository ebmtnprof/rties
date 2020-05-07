
######### The "rtiesShared" file includes functions that support all the rties analyses


################ Data manipulation functions


########### dataPrep

#' Reformat a user-provided dataframe in a generic form appropriate for \emph{rties} modeling
#'
#' The dataframe must be in a specific format and include several specific variables. See the "overview_data_prep" vignette for complete details on the necessary format and follow it closely if you'd like to avoid error messages. That vignette also includes information on how to structure the data if you have two variables within people (rather thean two people within dyads) or have indistinguishable dyads.
#'
#' @param basedata A user-provided dataframe that includes all variables needed for an rties analysis.
#' @param dyadId The name of the column in the dataframe that has the dyad-level identifier.
#' @param personId The name of the column in the dataframe that has the person-level identifier.
#' @param obs_name The name of the column in the dataframe that has the time-varying observable (e.g., the variable for which dynamics will be assessed).
#' @param dist_name The name of the column in the dataframe that has a variable that distinguishes the partners (e.g., sex, mother/daughter, etc) that is numeric and scored 0/1. 
#' @param time_name The name of the column in the dataframe that indicates sequential temporal observations.
#' @param time_lag An optional argument for the number of lags for the lagged observable. If a number is provided, the observed variable is lagged that amount. The other option is to use "absMaxCC". In this case the maximum cross-correlation is found for each dyad and the lag at which that occurs is used to lag their observed variables.
#'
#' @examples 
#' data <- rties_ExampleDataShort
#' newData <- dataPrep(basedata=data, dyadId="couple", personId="person", obs_name="dial", 
#' dist_name="female", time_name="time", time_lag=2)
#' head(newData) 
#'  
#' @return The function returns a dataframe that has all the variables needed for modeling system dynamics, each renamed to a generic variable name, which are:
#' \itemize{
#' \item id = person id
#' \item dyad = dyad id
#' \item obs = observed state variable
#' \item dist1 = 0/1 variable where the 1's indicate the 1's in the original distinguishing variable
#' \item time = the variable indicating temporal sequence
#' \item dist0 = 0/1 variable where the 1's indicate the 0's in the original distinguishing variable
#' \item obs_deTrend = the observed state variable with each person's linear trend removed
#' \item p_ = all the same variables, but for a person's partner rather than themselves
#'}

#' @export

dataPrep <- function(basedata, dyadId, personId, obs_name, dist_name, time_name, time_lag=NULL){
  
  vars <- c(dyadId, personId, obs_name, dist_name, time_name)
  basedata <- basedata[vars]
  names(basedata) <- c("dyad","id","obs","dist1","time")
    
      # check distinguishing variable is numeric 
  if (!is.numeric(basedata$dist1)){
	stop("the distinguishing variable must be a 0/1 numeric variable")
  }

  # create the dist0 variable
  basedata$dist0 <- ifelse(basedata$dist1 == 1, 0, 1)
    
  # check partners have same number of observations 
  notEqual <- vector()
  t <- table(basedata$dist1, basedata$dyad)
  for(i in 1:ncol(t)){		
    if (t[1,i] != t[2,i]){
	notEqual <- append(notEqual, as.numeric(dimnames(t)[[2]][i]))
	}
  }				
	if (length(notEqual) > 0){
	  print(notEqual)
	  stop("the partners in these dyads have unequal number of observations")
      notEqual <- NULL
	}
	    
  # center each person's data around their own regression line
  basedata <- lineCenterById(basedata)
  
  # create lagged variables
  if(!is.null(time_lag)){
	lag <- time_lag
	  
  if(lag == "absMaxCC"){
	crossCorr <- makeCrossCorr(basedata=basedata, dyadId="dyad", personId="id", obs_name="obs_deTrend", dist_name="dist1")
	cc <- crossCorr[!duplicated(crossCorr$dyad), ]
    lagTemp <- as.numeric(unlist(cc$maxLag))
    lag <- ifelse(lagTemp == 0, 1, lagTemp)
    dID <- unique(factor(basedata$dyad))
    lagData <- list()

      for (i in 1:length(dID)){
        lagi <- lag[i]
        datai <- basedata[basedata$dyad == dID[i], ]
        datai <- suppressMessages(DataCombine::slide(datai, Var="obs_deTrend", GroupVar="id", NewVar="obs_deTrend_Lag", slideBy= lagi))
        lagData[[i]] <- datai
       }
         basedata <- as.data.frame(do.call(rbind, lagData))
   } else {
	  basedata <- suppressMessages(DataCombine::slide(basedata, Var="obs_deTrend", GroupVar="id", NewVar="obs_deTrend_Lag", slideBy= -lag))
	 }
  }
  
  # put data in actor-partner format
  basedata <- actorPartnerDataTime(basedata, "dyad", "id")  
  return(basedata)
}


################## makeDist

#' Create a distinguishing variable (called "dist") for non-distinguishable dyads by assigning the partner who is lower on a chosen variable a 0 and the partner who is higher on the variable a 1. 
#'
#' @param basedata A user-provided dataframe.
#' @param dyadId The name of the column in the dataframe that has the couple-level identifier.
#' @param personId The name of the column in the dataframe that has the person-level identifier.
#' @param time_name The name of the column in the dataframe that indicates sequential temporal observations.
#' @param dist_name The name of the column in the dataframe that holds the variable to use for distinguishing the partners. For example, if "influence" was the variable, for each dyad the partner scoring lower on "influence" would be given a score of 0 on "dist" and the partner scoring higher on "influence" would be given a score of 1 on "dist"
#' @examples
#' data <- rties_ExampleDataShort
#' newData <- makeDist(basedata=data, dyadId="couple", personId="person", 
#' time_name="time", dist_name="relstress")
#' summary(newData$dist)


#' @return The function returns the original dataframe with an additional variable, called "dist" that distinguishes between partners based on the user-specified variable indicated by "dist_name"

#' @export

makeDist <- function(basedata, dyadId, personId, time_name, dist_name)
{
    vars1 <- c(dyadId, personId, time_name, dist_name)
    temp1 <- basedata[vars1]
    temp2 <- actorPartnerDataTime(temp1, dyadId, personId)     
    temp2$dist <- ifelse(temp2[ ,4] == temp2[ ,8], NA, 
				ifelse(temp2[ ,4] < temp2[ ,8], 0, 1))
    vars2 <- c("dyad", "person", "time", "dist")
    temp3 <- temp2[vars2]
    colnames(temp3)[colnames(temp3) == "dyad"] <- dyadId
    colnames(temp3)[colnames(temp3) == "person"] <- personId
    colnames(temp3)[colnames(temp3) == "time"] <- time_name
     
    temp4 <- plyr::join(basedata, temp3)
    return(temp4)
 }

############### lineCenterById

# This function creates a person-centered version of the observed variable "obs" called "obs_deTrend" that is centered around a linear regression line for each person, e.g., "obs_deTrend" is the residuals of "obs" when predicted from a linear regression on "time" for each person one at a time. 

lineCenterById <- function(basedata)
{
  newId <- unique(factor(basedata$id))
  dataCent <- list()
  for(i in 1:length(newId)){
	datai <- basedata[basedata$id == newId[i],]
	datai$obs_deTrend <- stats::resid(stats::lm(obs ~ time, data=datai, na.action=na.exclude))
	dataCent[[i]] <- datai
  }		
  basedata <- as.data.frame(do.call(rbind, dataCent)) 	
}		

########### removeDyads

#' Remove data for specified dyads from a dataframe
#'
#' Useful for cleaning data if some dyads have extensive missing or otherwise problematic data.
#'
#' @param basedata A user provided dataframe.
#' @param dyads A vector of dyad IDs to remove.
#' @param dyadId The variable in the dataframe specifying dyad ID.
#' @examples
#' data <- rties_ExampleDataShort
#' dyads <- c(3, 12)
#' newData <- removeDyads(basedata=data, dyads=dyads, dyadId="couple")
#'
#' @return A dataframe with the data for the specified dyads removed.

#' @export

removeDyads <- function (basedata, dyads, dyadId){
	colnames(basedata)[colnames(basedata)== dyadId] <- "dyad"
	basedata <- basedata[!basedata$dyad %in% dyads, ]
	colnames(basedata)[colnames(basedata)== "dyad"] <- dyadId
	return(basedata)
}

################# actorPartnerDataCross

#' Takes individual cross-sectional data from dyads and turns it into actor-partner format.
#'
#' Need to use a person ID that has first person in dyad numbered 1-n and second person in dyad = ID + some number larger than the number of dyads. Need dyad ID numbered same as for person ID for the first person in the dyad. Both members in each dyad need to have the same number of rows (rows of missing data are ok)
#'
#' @param basedata A dataframe with cross-sectional dyadic data.
#' @param dyadId The name of variable indicating dyad ID.
#' @param personId The name of the variable indicating peron ID.
#' @examples
#' data <- rties_ExampleDataShort
#' newData1 <- data[data$time==1, ] # make a cross-sectional dataframe
#' newData2 <- actorPartnerDataCross(basedata=newData1, dyadId="couple", personId="couple")
#' head(newData2)
#'
#' @return A dataframe in actor-partner format.

#' @export

actorPartnerDataCross <- function(basedata, dyadId, personId){
	
	colnames(basedata)[colnames(basedata) == dyadId] <- "dyad"
	colnames(basedata)[colnames(basedata) == personId] <- "person"	
	basedata <- basedata[order(basedata$person), ]
	dataA <- basedata
		
	P1 <- basedata[ which(basedata$dyad == basedata$person), ]
	P2 <- basedata[ which(basedata$dyad != basedata$person), ]
	
	P1_part <- P2
	P2_part <- P1
	colnames(P1_part) <- paste("p", colnames(P1_part), sep="_")
	colnames(P2_part) <- paste("p", colnames(P2_part), sep="_")
	dataP <- rbind(P1_part, P2_part)
	dataAP <- cbind(dataA, dataP)
	
	return(dataAP)		
}

################# actorPartnerDataTime

#' Takes individual repeated measures data from dyads and turns it into actor-partner format.
#'
#' Need to use a person ID that has first person in dyad numbered 1-n and second person in dyad = ID + some number larger than the number of dyads. Need dyad ID numbered same as for person ID for the first person in the dyad. Both members in each dyad need to have the same number of rows (rows of missing data are ok).
#'
#' @param basedata A dataframe with repeated measures dyadic data
#' @param dyadId The name of variable indicating dyad ID.
#' @param personId The name of the variable indicating peron ID.
#' @examples
#' data <- rties_ExampleDataShort
#' newData <- actorPartnerDataTime(basedata=data, dyadId="couple", personId="couple")
#' head(newData)
#'
#' @return A dataframe in actor-partner format.

#' @export

actorPartnerDataTime <- function(basedata, dyadId, personId){
		
    	colnames(basedata)[colnames(basedata) == dyadId] <- "dyad"
	colnames(basedata)[colnames(basedata) == personId] <- "person"
	basedata <- basedata[order(basedata$person), ]
	dataA <- basedata

    dID <- unique(factor(basedata$dyad))
	dataAP <- list()

	for(i in 1:length(dID))
		{
		datai <- basedata[basedata$dyad == dID[i],]
		dataA <- datai
		
		P1 <- datai[ which(datai$dyad == datai$person), ]
	    P2 <- datai[ which(datai$dyad != datai$person), ]

		P1_part <- P2
		P2_part <- P1
		colnames(P1_part) <- paste("p", colnames(P1_part), sep="_")
		colnames(P2_part) <- paste("p", colnames(P2_part), sep="_")
		dataP <- rbind(P1_part, P2_part)
		APi <- cbind(dataA, dataP)
		dataAP[[i]] <- APi		
		}		
	
	dataAP <- as.data.frame(do.call(rbind, dataAP))

} 	

############# makeFullData

#' Combines profile membership data from the latent profile analysis with other data for using the profile membership to predict and be predicted by the system variable.
#'
#' @param basedata The original dataframe provided by the user that includes all variables needed for an rties analysis, including potential system and control variables, etc.
#' @param dyadId The name of the column in the dataframe that has the couple-level identifier.
#' @param personId The name of the column in the dataframe that has the person-level identifier.
#' @param dist_name The name of the column in the dataframe that has a variable that distinguishes the partners (e.g., sex, mother/daughter, etc) that is numeric and scored 0/1. 
#' @param lpaData The object returned by the "inspectProfiles" function
#' @param params The list called "params" returned by one of the "indiv" functions (e.g., indivInertCoord or indivClo) 
#' @examples
#' data <- rties_ExampleDataShort
#' newData <- dataPrep(basedata=data, dyadId="couple", personId="person", 
#' obs_name="dial", dist_name="female", time_name="time", time_lag=2)
#' ic <- indivInertCoord(prepData=newData, whichModel="inertCoord")
#' profiles <- inspectProfiles(whichModel="inertCoord", prepData=newData, 
#' paramEst=ic$params, n_profiles=2)
#' fullData <- makeFullData(basedata=data, dyadId="couple", personId="person", 
#' dist_name="female", lpaData=profiles, params=ic$params)
#' head(fullData)
#'
#' @return A dataframe that contains all variables needed for using the profiles to predict, or be predicted by, the system variable.

#' @export

makeFullData <- function(basedata, dyadId, personId, dist_name, lpaData, params){
  colnames(basedata)[colnames(basedata)== dyadId] <- "dyad"
  colnames(basedata)[colnames(basedata)== personId] <- "person"
  colnames(basedata)[colnames(basedata)== dist_name ] <- "dist1"

  basedata <- basedata[!duplicated(basedata$person), ]
  basedata$dist0 <- ifelse(basedata$dist1 == 1, 0, 1)

  fullData <- plyr::join(basedata, lpaData)
  
  return(fullData)
}

############# Max_Min_CCF_Signed

#' A helper function for makeCrossCorr

#' @param a First time-series used in the cross-correlation
#' @param b Second time-series used in the cross-correlation
#' 
#' @return A list of maximum absolute value cross-correlations and the lag at which they occurred.

#' @importFrom stats ccf na.exclude

Max_Min_CCF_Signed <- function (a, b) {
d <- ccf(a, b, plot = FALSE, na.action=na.exclude)
cor = d$acf[ ,,1] 
lag = d$lag[ ,,1] 
res = data.frame(cor,lag) 
for(i in 1:nrow(res)) if(is.na(res[i,1])) max=NA else {
res_max = res[which.max(res$cor),] 
res_min = res[which.min(res$cor),]
res_min_abs = abs(res_min)
res.max.abs = which.max(c(res_max$cor, res_min_abs$cor))
if (res.max.abs==1) max=res_max else max=res_min}
output <-data.frame(max = max)
output
}


############# makeCrossCorr

#' Calculates cross-correlations for a given variable and returns a dataframe with the largest absolute cross-correlation and its lag added for each dyad (e.g., it returns either the most negative or most positive cross-correlation, whichever is larger in absolute terms).
#'
#' @param basedata The original dataframe provided by the user that includes all variables needed for an rties analysis, including potential system and control variables, etc.
#' @param dyadId The name of the column in the dataframe that has the couple-level identifier.
#' @param personId The name of the column in the dataframe that has the person-level identifier.
#' @param obs_name The name of the column in the dataframe that has the time-varying observable (e.g., the variable for which dynamics will be assessed).
#' @param dist_name The name of the column in the dataframe that has a variable that distinguishes the partners (e.g., sex, mother/daughter, etc) that is numeric and scored 0/1. 
#' @examples
#' data <- rties_ExampleDataShort
#' newData <- makeCrossCorr(basedata=data, dyadId="couple", personId="person", 
#' obs_name="dial", dist_name="female")
#' head(newData)
#'
#' @return A cross-sectional version of the original dataframe with maximal absolute-value cross-correlations and their lags added for each dyad.

#' @export

makeCrossCorr <- function(basedata, dyadId, personId, obs_name, dist_name){
  
  newdata <- basedata
  colnames(newdata)[colnames(newdata)== dyadId] <- "dyad"
  colnames(newdata)[colnames(newdata)== personId] <- "person"
  colnames(newdata)[colnames(newdata)== obs_name] <- "dv"
  colnames(newdata)[colnames(newdata)== dist_name ] <- "dist1"

  crossCorr <- list()
  dID <- unique(factor(newdata$dyad))

  for (i in 1:length(dID)){
  datai <- newdata[newdata$dyad == dID[i], ]
  dist1 <- datai[ which(datai$dist1 == 1), "dv"]
  dist0 <- datai[ which(datai$dist1 == 0), "dv"]
  crossCorr[[i]] <- Max_Min_CCF_Signed(dist0, dist1)
  }

  cc1 <- lapply(crossCorr, function(x) lapply(x, function(x) ifelse(is.null(x), NA, x)))
  cc2 <- lapply(cc1, function(x) lapply(x, function(x) ifelse(is.numeric(x), round(x,digits=3), x)))
  cc <- as.data.frame(do.call(rbind, crossCorr))
  cc$newID <- dID
  colnames(cc) <- c("maxCor","maxLag", dyadId)
 
  temp1 <- newdata[!duplicated(newdata$person), ]
  colnames(temp1)[colnames(newdata)== "dyad"] <- dyadId
  colnames(temp1)[colnames(newdata)== "person"] <- personId
  temp2 <- plyr::join(cc, temp1)
  
  crossCorr <- temp2
  return(crossCorr)
 }


################ Plotting functions

########## histAll

#' Histograms for all numeric variables in a dataframe.
#'
#' Useful for checking distributions to assess normality
#'
#' @param basedata A user-provided dataframe.
#' @examples
#' data <- rties_ExampleDataShort
#' vars <- c("reltime","ambiv","love","conflict")
#' newData <- data[vars ]
#' histAll(newData)
#' 
#' @return No return value. Prints plots to the console.

#' @export

histAll <- function(basedata)
{
  nums <- sapply(basedata, is.numeric)
  numdata <- basedata[ ,nums]
  
  opar <- graphics::par(no.readonly =TRUE) 
  on.exit(graphics::par(opar))

  graphics::par(mfrow=c(2,2))
  for(i in 1:length(numdata)){
	graphics::hist(numdata[,i], main=NULL, xlab=names(numdata[i]))
  }
}

############### plotRaw

#' Plots of observed variable over time by dyad.
#'
#' Produces plots of the observed variable for each dyad over time to check for data errors, etc. 
#'
#' @param basedata A user provided dataframe.
#' @param dyadId The name of the column in the dataframe that has the dyad-level identifier.
#' @param obs_name The name of the column in the dataframe that has the time-varying observable (e.g., the variable for which dynamics will be assessed).
#' @param dist_name The name of the column in the dataframe that has a variable that distinguishes the partners (e.g., sex, mother/daughter, etc) that is numeric and scored 0/1.
#' @param time_name The name of the column in the dataframe that indicates sequential temporal observations.
#' @param dist0name An optional name for the level-0 of the distinguishing variable to appear on plots (e.g., "Women").
#' @param dist1name An optional name for the level-1 of the distinguishing variable to appear on plots (e.g., "Men").
#' @param plot_obs_name An optional name for the observed state variable to appear on plots (e.g., "Emotional Experience").
#' @param printPlots If true (the default) plots are displayed on the screen.
#' @examples
#' data <- rties_ExampleDataShort
#' plotRaw(basedata=data, dyad="couple", obs_name="dial", dist_name="female", time_name="time")

#' @return A list of plots.

#' @export

plotRaw <- function(basedata, dyadId, obs_name, dist_name, time_name, dist0name=NULL, dist1name=NULL, plot_obs_name=NULL, printPlots=T) 
{
  basedata <- basedata[ ,c(dyadId, obs_name, dist_name, time_name) ]
  names(basedata) <- c("dyad", "obs", "dist", "time")
  
  # check distinguishing variable is numeric 
  if (!is.numeric(basedata$dist)){
	stop("the distinguishing variable must be a 0/1 numeric variable")
  }

  if(is.null(dist0name)){dist0name <- "dist0"}
  if(is.null(dist1name)){dist1name <- "dist1"}
  if(is.null(plot_obs_name)){plot_obs_name <- "obs"}
  
  dist <- NULL
  plots <- lattice::xyplot(obs~time|as.factor(dyad), data = basedata, group=dist, type=c("l"), ylab=plot_obs_name, col=c("red", "blue"), key=list(space="right", text=list(c(dist1name,dist0name)), col=c("blue", "red")),as.table=T, layout = c(3,3))
  
  if(printPlots==T){print(plots)}
  
  return(plots)
}

################ plotDataByProfile

#' Plots of de-trended observed variable over time, with dyads separated into groups based on LPA profile membership.
#'
#' @param prepData A dataframe created by the dataPrep function.
#' @param fullData A dataframe created by the makeFullData function.
#' @param n_profiles The number of profiles that were estimated.
#' @param dist0name An optional name for the level-0 of the distinguishing variable (e.g., "Women"). Default is dist0.
#' @param dist1name An optional name for the level-1 of the distinguishing variable (e.g., "Men"). Default is dist1.
#' @param plot_obs_name An optional name for the observed state variable to appear on plots (e.g., "Emotional Experience").
#' @param printPlots If true (the default) plots are displayed on the screen.
#' @examples
#' # See vignettes for examples.

#' @return A list of plots.

#' @export

plotDataByProfile <- function(prepData, fullData, n_profiles, dist0name=NULL, dist1name=NULL, plot_obs_name=NULL, printPlots=T){

  if(is.null(dist0name)){dist0name <- "dist0"}
  if(is.null(dist1name)){dist1name <- "dist1"}
  if(is.null(plot_obs_name)){plot_obs_name <- "obs_deTrend"}
    
    vars1 <- c("dyad", "dist0","profile")
    temp1 <- fullData[ , vars1]
    vars2 <- c("dyad", "dist0", "obs_deTrend", "time")
    temp2 <- prepData[ , vars2]
    
    dist0 <- NULL
    temp3 <- plyr::join(temp1, temp2)

    plots <- list()
    for(i in 1:n_profiles){
      tempi <- temp3[ which(temp3$profile == i), ]
      label <- paste("Profile", i, sep="-")   
      plots[[i]] <- lattice::xyplot(obs_deTrend ~ time|as.factor(dyad), data = tempi, group=dist0, type=c("l"), ylab=plot_obs_name, main=label, col=c("red", "blue"), key=list(space="right", text=list(c(dist1name,dist0name)), col=c("blue", "red")), as.table=T, layout = c(3,3))
      }
      
    if(printPlots==T){print(plots)}
    
    return(plots)
}

