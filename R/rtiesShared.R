

######### The "rtiesShared" file includes functions that support all the rties analyses

#' Reformat a user-provided dataframe in a generic form appropriate for \emph{rties} modeling
#'
#' In the dataframe, the partners within each dyad must have the same number of observations (e.g. rows of data), although those can include rows that have missing values (NAs). Each dyad, however, can have it's own unique number of observations.
#'
#' @param basedata A user-provided dataframe.
#' @param id The name of the column in the dataframe that has the person-level identifier.
#' @param dyad The name of the column in the dataframe that has the dyad-level identifier.
#' @param obs The name of the column in the dataframe that has the time-varying observable (e.g., the variable for which dynamics will be assessed).
#' @param sysVar The name of the column in the dataframe that has the system variable (e.g., something that will be predicted from the dynamics of the system).
#' @param dist The name of the column in the dataframe that has a variable that distinguishes the partners (e.g., sex, mother/daughter, etc) that is numeric and scored 0/1.
#' @param time_name The name of the column in the dataframe that indicates sequential temporal observations.
#' @param time_lag An optional argument for the number of lags for the lagged observable.
#' @param center. An optional vector of centering values for the system variable, for use with versions of the models that use it as a moderator, rather than an outcome.
#'
#' @return The function returns a dataframe that has all the variables needed for rties modeling, each renamed to a generic variable name, which are:
#' \itemize{
#' \item id = person id
#' \item dyad = dyad id
#' \item obs = observed state variable
#' \item sysVar = system variable
#' \item dist1 = 0/1 variable where the 1's indicate the 1's in the original distinguishing variable
#' \item time = the variable indicating temporal sequence
#' \item dist0 = 0/1 variable where the 1's indicate the 0's in the original distinguishing variable
#' \item obs_deTrend = the observed state variable with each person's linear trend removed
#' \item p_ = all the same variables, but for a person's partner rather than themselves
#' \item if the centering option is used, it will also return sysVarL, sysVarM, and sysVarH which are centered versions of the system variable}

#' @export
dataPrep <- function(basedata,id,dyad,obs,sysVar,dist,time_name,time_lag=NULL, center=NULL) 
{
  basedata <- subset(basedata, select=c(id, dyad, obs, sysVar, dist, time_name))
  names(basedata)[1] <- "id"
  names(basedata)[2] <- "dyad"
  names(basedata)[3] <- "obs"
  names(basedata)[4] <- "sysVar"
  names(basedata)[5] <- "dist1"
  names(basedata)[6] <- "time"
  
     # check distinguishing variable is numeric 
    if (!is.numeric(basedata$dist1))		
	  {
		cat("\n error: the distinguishing variable must be a 0/1 numeric variable\n")
		stop(call.=F)
    }

    ## create the dist0 variable
    basedata$dist0 <- ifelse(basedata$dist1 == 1, 0, 1)
    
   # check partners have same number of observations and exit with error message if not
    	notEqual <- vector()
    	t <- table(basedata$dist1, basedata$dyad)
    	for(i in 1:ncol(t))
    	{		if (t[1,i] != t[2,i])
			{
			notEqual <- append(notEqual, as.numeric(dimnames(t)[[2]][i]))
			}
		}				
	if (length(notEqual) > 0)		
	{
		cat("\n error: the partners in these dyads have unequal number of observations\n")
		print(notEqual)
		stop(call.=F)
		rm(notEqual, envir = .GlobalEnv)
	}
	    
	 basedata <- lineCenterById(basedata)

	 if(!is.null(center))
	 {
	    basedata$sysVarL <- basedata$sysVar - center[1]
	    basedata$sysVarM <- basedata$sysVar - center[2]
	    basedata$sysVarH <- basedata$sysVar - center[3]
	 }
	   
	  if(!is.null(time_lag))
	  {
	   lag <- time_lag
	   basedata <- DataCombine::slide(basedata, Var="obs_deTrend", GroupVar="id", NewVar="obs_deTrend_Lag", slideBy= -lag)
	  }
	    basedata <- actorPartnerDataTime(basedata, basedata$dyad, basedata$id)
	    
	   return(basedata)
}

############### lineCenterById

# This function creates a person-centered version of the observed variable "obs" called "obs_deTrend" that is centered around a linear regression line for each person, e.g., "obs_deTrend" is the residuals of "obs" when predicted from a linear regression on "time" for each person one at a time. 

lineCenterById <- function(basedata)
{
	newId <- unique(factor(basedata$id))
	dataCent <- list()
	for(i in 1:length(newId))
	{
		datai <- basedata[basedata$id == newId[i],]
		datai$obs_deTrend <- resid(lm(obs ~ time, data=datai, na.action=na.exclude))
		dataCent[[i]] <- datai
	}		
	basedata <- as.data.frame(do.call(rbind, dataCent)) 	
}		


############# plotting functions: The following are a set of useful basic plots

#' Histograms for all numeric variables in a dataframe.
#'
#' Useful for checking distributions of potential system variables to assess normality
#'
#' @param basedata A user-provided dataframe.

#' @export
histAll <- function(basedata)
{
	nums <- sapply(basedata, is.numeric)
	numdata <- basedata[ ,nums]
	par(mfrow=c(4,4))
	for(i in 1:length(numdata))
	{
		hist(numdata[,i], xlab=NULL, main=names(numdata[i]))
	}
}

#' Plots of observed variable over time by dyad.
#'
#' Produces plots of the observed variable for each dyad over time to check for data errors, etc. 
#'
#' @param basedata A dataframe.
#' @param dyad The name of the column in the dataframe that has the dyad-level identifier.
#' @param obs The name of the column in the dataframe that has the time-varying observable (e.g., the variable for which dynamics will be assessed).
#' @param dist The name of the column in the dataframe that has a variable that distinguishes the partners (e.g., sex, mother/daughter, etc) that is numeric and scored 0/1.
#' @param time_name The name of the column in the dataframe that indicates sequential temporal observations.
#' @param dist0name A name for the 0-level of the distinguishing variable (e.g., "Women").
#' @param dist1name A name for the 1-level of the distinguishing variable (e.g., "Men").

#' @export

plotRaw <- function(basedata,dyad,obs,dist,time_name, dist0name, dist1name) 
{
  obs_name <- obs 
  basedata <- subset(basedata, select=c(dyad, obs, dist, time_name))
  names(basedata)[1] <- "dyad"
  names(basedata)[2] <- "obs"
  names(basedata)[3] <- "dist"
  names(basedata)[4] <- "time"
  
 lattice::xyplot(obs~time|as.factor(dyad), data = basedata, group=dist, type=c("l"), ylab=obs_name, col=c("red", "blue"), key=list(space="right", text=list(c(dist1name,dist0name)), col=c("blue", "red")),as.table=T, layout = c(5,5))
  }

# Plots of linear regression lines for both people in each dyad
#'
#' Produces plots of temporal trajectories predicted by linear dyadic growth models. 
#'
#' @param basedata A dataframe produced by "dataPrep".
#' @param dist0name A name for the 0-level of the distinguishing variable (e.g., "Women").
#' @param dist1name A name for the 1-level of the distinguishing variable (e.g., "Men").
#' @param obsName A name for the observed variable being plotted (e.g., "Emotional Experience").

#' @export
plotLinear <- function(basedata, dist0name, dist1name, obsName){
	lattice::xyplot(obs~time|as.factor(dyad), data = basedata, group=dist1, type=c("r"), ylab=obsName, col=c("red", "blue"), key=list(space="right", text=list(c(dist1name,dist0name)), col=c("blue", "red")),as.table=T, layout = c(5,5))
}

# This function plots curvlinear (loess smoothed) lines for each person in each dyad

plotCurve <- function(basedata, dist0name, dist1name, obsName){
	lattice::xyplot(obs~time|as.factor(dyad), data = basedata, group=dist1, type=c("smooth"), ylab=obsName, col=c("red", "blue"), key=list(space="right", text=list(c(dist1name,dist0name)), col=c("blue", "red")),as.table=T, layout = c(5,5))
}

### orderedPlotsLinearAve: This function produces plots of the linear fits of the observed variable for each dyad in ascending order of the dyad averages on the sysVar variable, which is jittered to deal with matching values. 

orderedPlotsLinearAve <- function(basedata)
{
	temp <- aggregate(basedata$sysVar, by=list(basedata$dyad), FUN="mean", na.rm=TRUE)
	colnames(temp) <- c("dyad", "sysVarAve")
	temp$sysVarAveJ <- round(jitter(temp$sysVarAve), 5)
	newData <- join(basedata, temp)

	lattice::xyplot(obs~time|as.factor(sysVarAveJ), data = newData,group=dist1,type=c("r"), col=c("blue", "red"),
			as.table=T, layout = c(3,3))
}


## orderedPlotsLinearDist: This function produces plots of the linear fits of the observed variable for each dyad in ascending order of one of the partner's scores on the sysVar variable, which is jittered to deal with matching values. Which partner is used is determined by the dist argument, which must be zero or one to correspond with the options for the distinguishing variable

orderedPlotsLinearDist <- function(basedata, dist)
{
	if(dist==0)
	{
	temp <- basedata[basedata$time == 1 & basedata$dist0 == 1, ]
	} else if(dist==1)
			{temp <- basedata[basedata$time == 1 & basedata$dist1 == 1, ]}
			 else {cat("\n error: dist must be 0 or 1\n")
			 	stop(call.=F)}
	temp$sysVarJ <- round(jitter(temp$sysVar), 5)
	temp <- subset(temp, select=c(dyad, sysVarJ))
	data <- plyr::join(basedata, temp)
	
	lattice::xyplot(obs ~time|as.factor(sysVarJ), data = data, group=dist1, type=c("r"), col=c("red", "blue"),
			as.table=T, layout = c(3,3))
}


### orderedPlotsDetrendAve: This function produces plots of the trajectories of the detrended observed variable for each dyad in ascending order of the dyad averages on the sysVar variable, which is jittered to deal with matching values. 

orderedPlotsDetrendAve <- function(basedata)
{
	temp <- aggregate(basedata$sysVar, by=list(basedata$dyad), FUN="mean", na.rm=TRUE)
	colnames(temp) <- c("dyad", "sysVarAve")
	temp$sysVarAveJ <- round(jitter(temp$sysVarAve), 5)
	data <- plyr::join(basedata, temp)

	lattice::xyplot(obs_deTrend ~time|as.factor(sysVarAveJ), data = data, group=dist1, type=c("l"), col=c("red", "blue"),
			as.table=T, layout = c(3,3))
}


### orderedPlotsDetrendDist: This function produces plots of the trajectories of the detrended observed variable for each dyad in ascending order of one of the partner's scores on the sysVar variable, which is jittered to deal with matching values. Which partner is used is determined by the dist argument, which must be zero or one to correspond with the options for the distinguishing variable

orderedPlotsDetrendDist <- function(basedata, dist)
{
	if(dist==0)
	{
	temp <- basedata[basedata$time == 1 & basedata$dist0 == 1, ]
	} else if(dist==1)
			{temp <- basedata[basedata$time == 1 & basedata$dist1 == 1, ]}
			 else {cat("\n error: dist must be 0 or 1\n")
			 	stop(call.=F)}
	temp$sysVarJ <- round(jitter(temp$sysVar), 5)
	temp <- subset(temp, select=c(dyad, sysVarJ))
	data <- plyr::join(basedata, temp)
	
	lattice::xyplot(obs_deTrend ~time|as.factor(sysVarJ), data = data, group=dist1, type=c("l"), col=c("red", "blue"),
			as.table=T, layout = c(3,3))
}

 

sysVarByParam <- function(paramData, colToPlot, sysVarName)
{
	ymin <- min(paramData$sysVar, na.rm=T)
	ymax <- max(paramData$sysVar, na.rm=T)
	par(mfrow=c(4,4))
	for(i in colToPlot)
	{
		xmin <- min(paramData[i], na.rm=T)
		xmax <- max(paramData[i], na.rm=T)
		plot(paramData$sysVar ~ paramData[,i], xlab="", ylab=sysVarName, ylim=c(ymin, ymax), xlim=c(xmin, xmax), main=names(paramData[i]))
	}	
}

########### removeDyads

#' Remove data for specified dyads from a dataframe
#'
#' Useful for cleaning data if some dyads have extensive missing or otherwise problematic data.
#'
#' @param basedata A dataframe produced by "dataPrep"
#' @param dyads A vector of dyad IDs to remove.
#' @param dyadID The variable in the dataframe specifying dyad ID; should be in the form dataframe_name$variable_name (e.g., data$couple).
#'
#' @return A dataframe with the data for the specified dyads removed.

#' @export
removeDyads <- function (basedata, dyads, dyadID){
	basedata <- subset(basedata, !dyadID %in% dyads)
	return(basedata)
}


## actorPartnerDataCross: This function takes individual cross-sectional data from dyads and turns it into actor-partner format

# Need to use a person ID that has first person in dyad numbered 1-n and second person in dyad = ID + some number larger than the number of dyads
# Need dyad ID numbered same as for person ID for the first person in the dyad
# Both members in each dyad need to have the same number of rows (rows of missing data are ok)

#  arguments: basedata = name of an R data frame containing original data, dyadID = name of variable indicating dyad ID, and personID = name of the variable indicating peron ID.
# function will return a data file in actor-partner format

# Example:
# dataAP <- actorPartnerDataCross(data, data$dyad, data$person)

actorPartnerDataCross <- function(basedata, dyadID, personID){
	basedata$d <- dyadID
	basedata$p <- personID
	dataA <- basedata
	
	P1 <- subset(basedata, basedata$d == basedata$p)
	P2 <- subset(basedata, basedata$d != basedata$p)
	P1_part <- P2
	P2_part <- P1
	colnames(P1_part) <- paste("p", colnames(P1_part), sep="_")
	colnames(P2_part) <- paste("p", colnames(P2_part), sep="_")
	dataP <- rbind(P1_part, P2_part)
	dataAP <- cbind(dataA, dataP)
	dataAP <- subset(dataAP, select=-c(d, p, p_d, p_p))
	return(dataAP)		
}


## actorPartnerDataTime: This function takes individual repeated measures data from dyads and turns it into actor-partner format

# Need to use a person ID that has first person in dyad numbered 1-n and second person in dyad = ID + some number larger than the number of dyads
# Need dyad ID numbered same as for person ID for the first person in the dyad
# Both members in each dyad need to have the same number of rows (rows of missing data are ok)
#  arguments: basedata = name of an R data frame containing original data, dyadID = name of variable indicating dyad ID, and personID = name of the variable indicating peron ID.
# function will create a data file in actor-partner format

# Example:
# dataAP <- actorPartnerDataTime(data, data$Dyad, data$Person)

actorPartnerDataTime <- function(basedata, dyadID, personID)
{
    basedata$d <- dyadID
    basedata$p <- personID
    dID <- unique(factor(dyadID))
	dataAP <- list()

	for(i in 1:length(dID))
		{
		datai <- basedata[basedata$d == dID[i],]
		dataA <- datai
		P1 <- subset(datai, datai$d == datai$p)
		P2 <- subset(datai, datai$d != datai$p)
		P1_part <- P2
		P2_part <- P1
		colnames(P1_part) <- paste("p", colnames(P1_part), sep="_")
		colnames(P2_part) <- paste("p", colnames(P2_part), sep="_")
		dataP <- rbind(P1_part, P2_part)
		APi <- cbind(dataA, dataP)
		APi <- subset(APi, select=-c(d, p, p_d, p_p))
		dataAP[[i]] <- APi		
		}		
	dataAP <- as.data.frame(do.call(rbind, dataAP))
} 	

