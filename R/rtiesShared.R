

######### The "rtiesShared" file includes functions that support all the rties analyses


################ Data manipulation functions

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
#' @param center An optional vector of centering values for the system variable.
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
#'}

#' @export
dataPrep <- function(basedata,id,dyad,obs,sysVar,dist,time_name,time_lag=NULL, center=NULL){
  
  basedata <- subset(basedata, select=c(id, dyad, obs, sysVar, dist, time_name))
  names(basedata)[1] <- "id"
  names(basedata)[2] <- "dyad"
  names(basedata)[3] <- "obs"
  names(basedata)[4] <- "sysVar"
  names(basedata)[5] <- "dist1"
  names(basedata)[6] <- "time"
  
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
      rm(notEqual, envir = .GlobalEnv)
	}
	    
  basedata <- lineCenterById(basedata)

  if(!is.null(center)){
	basedata$sysVarL <- basedata$sysVar - center[1]
	basedata$sysVarM <- basedata$sysVar - center[2]
	basedata$sysVarH <- basedata$sysVar - center[3]
  }
	   
  if(!is.null(time_lag)){
	lag <- time_lag
	basedata <- suppressMessages(DataCombine::slide(basedata, Var="obs_deTrend", GroupVar="id", NewVar="obs_deTrend_Lag", slideBy= -lag))
  }
  
  basedata <- actorPartnerDataTime(basedata, "dyad", "id")  
  return(basedata)
}

############### lineCenterById

# This function creates a person-centered version of the observed variable "obs" called "obs_deTrend" that is centered around a linear regression line for each person, e.g., "obs_deTrend" is the residuals of "obs" when predicted from a linear regression on "time" for each person one at a time. 

lineCenterById <- function(basedata)
{
  newId <- unique(factor(basedata$id))
  dataCent <- list()
  for(i in 1:length(newId)){
	datai <- basedata[basedata$id == newId[i],]
	datai$obs_deTrend <- resid(lm(obs ~ time, data=datai, na.action=na.exclude))
	dataCent[[i]] <- datai
  }		
  basedata <- as.data.frame(do.call(rbind, dataCent)) 	
}		

########### removeDyads

#' Remove data for specified dyads from a dataframe
#'
#' Useful for cleaning data if some dyads have extensive missing or otherwise problematic data.
#'
#' @param basedata A dataframe.
#' @param dyads A vector of dyad IDs to remove.
#' @param dyadID The variable in the dataframe specifying dyad ID; should be in the form dataframe_name$variable_name (e.g., data$couple).
#'
#' @return A dataframe with the data for the specified dyads removed.

#' @export
removeDyads <- function (basedata, dyads, dyadID){
	basedata <- subset(basedata, !dyadID %in% dyads)
	return(basedata)
}

#' actorPartnerDataCross: This function takes individual cross-sectional data from dyads and turns it into actor-partner format.
#'
#' Need to use a person ID that has first person in dyad numbered 1-n and second person in dyad = ID + some number larger than the number of dyads. Need dyad ID numbered same as for person ID for the first person in the dyad. Both members in each dyad need to have the same number of rows (rows of missing data are ok)
#'
#' @param basedata A dataframe with cross-sectional dyadic data.
#' @param dyadID The name of variable indicating dyad ID.
#' @param personID The name of the variable indicating peron ID.
#' @return A dataframe in actor-partner format.

# @export

actorPartnerDataCross <- function(basedata, dyadID, personID){
	
	basedata$d <- basedata[, dyadID]
	basedata$p <- basedata[, personID]
	basedata <- basedata[order(basedata$p), ]
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

#' actorPartnerDataTime: This function takes individual repeated measures data from dyads and turns it into actor-partner format.
#'
#' Need to use a person ID that has first person in dyad numbered 1-n and second person in dyad = ID + some number larger than the number of dyads. Need dyad ID numbered same as for person ID for the first person in the dyad. Both members in each dyad need to have the same number of rows (rows of missing data are ok).
#'
#' @param basedata A dataframe with repeated measures dyadic data
#' @param dyadID The name of variable indicating dyad ID.
#' @param personID The name of the variable indicating peron ID.
#' @return A dataframe in actor-partner format.

#' @export

actorPartnerDataTime <- function(basedata, dyadID, personID){
		
    basedata$d <- basedata[, dyadID]
    basedata$p <- basedata[, personID]
    basedata <- basedata[order(basedata$p), ]
    dID <- unique(factor(basedata$d))
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

#' Combines profile membership data from the latent profile analysis with other data for using the profile membership to predict and be predicted by the system variable.
#'
#' @param modelData A dataframe created by one of the "indiv" functions containing the parameter estimates for one of the models (e.g., inertCoord or clo) in combination and all other data in the analysis (e.g., sysVar, covariates, etc)
#' @param lpaData The object created by tidyLPA's "estimate_profiles" function when "return_orig_df = TRUE"
#' @param lpaParams The object created by tidyLPA's "estimate_profiles" function when "to_return = mclust"
#' @param whichModel The name of the model that is being investigated (e.g., "inertCoord" or "clo")
#' @return A list containing 1) a dataframe that contains all variables needed for plotting profile trajectories and using the profiles to predict, or be predicted by, the system variable (called "profileData"), and 2)a dataframe containing the profile parameter estimates needed for plotting the predicted trajectories for each profile (called "profileParams").

#' @export

makeLpaData <- function(modelData, lpaData, lpaParams, whichModel){
  
  if(whichModel !="inertCoord" & whichModel != "clo" ){
	stop("the model type must be either inertCoord or clo")
	
	} else if (whichModel == "inertCoord"){
		
	    data <- as.data.frame(lpaData)
        data <- subset(data, select=c(dyad, profile))
        data <- suppressMessages(plyr::join(modelData, data))
        data$profileN <- as.numeric(data$profile) - 1        
        params <- lpaParams$parameters$mean
	
	} else if (whichModel == "clo"){
			
	    data <- as.data.frame(lpaData)
        data <- subset(data, select=c(dyad, profile))
        data <- suppressMessages(plyr::join(modelData, data))
        data <- data[!duplicated(data$id), ] 
        data$profileN <- as.numeric(data$profile) - 1        
        params <- lpaParams$parameters$mean
	}
	
  results <- list(profileData=data, profileParams=params)
}




################ Plotting functions

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

  par(mfrow=c(2,2))
  for(i in 1:length(numdata)){
	hist(numdata[,i], main=NULL, xlab=names(numdata[i]))
  }
  par(mfrow=c(1,1))
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
#' @param dist0name An optional name for the level-0 of the distinguishing variable (e.g., "Women"). Default is dist0.
#' @param dist1name An optional name for the level-1 of the distinguishing variable (e.g., "Men"). Default is dist1.

#' @export

plotRaw <- function(basedata, dyad, obs, dist, time_name, dist0name=NULL, dist1name=NULL) 
{
  obs_name <- obs 
  basedata <- basedata[ ,c(dyad, obs, dist, time_name) ]
  names(basedata) <- c("dyad", "obs", "dist", "time")
  
  # check distinguishing variable is numeric 
  if (!is.numeric(basedata$dist)){
	stop("the distinguishing variable must be a 0/1 numeric variable")
  }

  if(is.null(dist0name)){dist0name <- "dist0"}
  if(is.null(dist1name)){dist1name <- "dist1"}
 
  lattice::xyplot(obs~time|as.factor(dyad), data = basedata, group=dist, type=c("l"), ylab=obs_name, col=c("red", "blue"), key=list(space="right", text=list(c(dist1name,dist0name)), col=c("blue", "red")),as.table=T, layout = c(3,3))
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

######################### Miscellaneous functions


biserialCor <- function (x, y, level = 1) 
{
    if (!is.numeric(x)) 
        stop("'x' must be a numeric variable.\n")
    y <- as.factor(y)
    if (length(levs <- levels(y)) > 2) 
        stop("'y' must be a dichotomous variable.\n")
    if (length(x) != length(y)) 
        stop("'x' and 'y' do not have the same length")
    
        cc.ind <- complete.cases(x, y)
        x <- x[cc.ind]
        y <- y[cc.ind]
    
    ind <- y == levs[level]
    diff.mu <- mean(x[ind]) - mean(x[!ind])
    prob <- mean(ind)
    sd.pop <- sd(x) * sqrt((length(x) - 1)/length(x))
    diff.mu * sqrt(prob * (1 - prob))/sd.pop
}
