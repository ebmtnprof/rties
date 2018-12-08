

######## This file includes all the functions needed for an inertia-coordination analysis


#' Produces auto-correlation plots of the observed state variable for lags of -+ 20 time steps for each dyad.
#' 
#' @param basedata A dataframe that has the variables needed for the dataPrep function, either in the right order or as named arguments (e.g., the order and names are: id, dyad, obs, sysVar, dist, time_name). See documentation for the dataPrep function for more detail.

#' @export

autoCorPlots <- function(basedata)
{
  names(basedata)[1] <- "id"
  names(basedata)[2] <- "dyad"
  names(basedata)[3] <- "obs"
  names(basedata)[4] <- "sysVar"
  names(basedata)[5] <- "dist1"
  names(basedata)[6] <- "time"

  basedata <- basedata[complete.cases(basedata), ]
  basedata <- lineCenterById(basedata)
  basedata <- actorPartnerDataTime(basedata, "dyad", "id")
  
  newDiD <- unique(factor(basedata$dyad))
  acf <- list()
  plotTitle <- list()
	
  for (i in 1:length(newDiD)){
	datai <- basedata[basedata$dyad == newDiD[i], ]
	datai_ts <- zoo::zoo(datai[,7])
	d <- acf(datai_ts, plot=F, na.action=na.exclude, lag.max=20)
    acf[[i]] <- d
    plotTitle[[i]] <- as.character(unique(datai$dyad))
  }  
  
  par(mfrow=c(3,3))
  for(j in 1:length(newDiD))
  plot(acf[[j]], main=paste("AutoCorr_Dyad", plotTitle[j], sep="_"))
  par(mfrow=c(1,1))
  return(acf)
}

#' Produces cross-correlation plots of the observed state variable for lags of -+ 20 time steps for each dyad.
#' 
#' @param basedata A dataframe that has the variables needed for the dataPrep function, either in the right order or as named arguments (e.g., the order and names are: id, dyad, obs, sysVar, dist, time_name). See documentation for the dataPrep function for more detail.

#' @export

crossCorPlots <- function(basedata)
{
  names(basedata)[1] <- "id"
  names(basedata)[2] <- "dyad"
  names(basedata)[3] <- "obs"
  names(basedata)[4] <- "sysVar"
  names(basedata)[5] <- "dist1"
  names(basedata)[6] <- "time"

  basedata <- basedata[complete.cases(basedata), ]
  basedata <- lineCenterById(basedata)
  basedata <- actorPartnerDataTime(basedata, "dyad", "id")
  
  newDiD <- unique(factor(basedata$dyad))
  ccf <- list()
  plotTitle <- list()
	
  for (i in 1:length(newDiD)){
	datai <- basedata[basedata$dyad == newDiD[i], ]
	datai_ts1 <- zoo::zoo(datai[ ,7])
	datai_ts2 <- zoo::zoo(datai[ ,14])
	d <- ccf(datai_ts1, datai_ts2, type="correlation", plot=F, na.action=na.exclude, lag.max=20)
    ccf[[i]] <- d
    plotTitle[[i]] <- as.character(unique(datai$dyad))
  }  
  
  par(mfrow=c(3,3))
  for(j in 1:length(newDiD))
  plot(ccf[[j]], main=paste("CrossCorr_Dyad", plotTitle[j], sep="_"))
  par(mfrow=c(1,1))
  return(acf)
}



#' Estimates versions of the inertia-coordination model for each dyad.
#' 
#' The user specifies which of 3 models are to be estimated. Each model predicts the observed state variables (with linear trends removed) from either: 1) Inertia only ("inert")- each person's intercept and each person's own observed state variable lagged at the amount specified during the dataPrep step (again with linear trends removed), 2) Coordination only ("coord")- each person's intercept and each person's partner's state variable lagged at the amount specified (again with linear trends removed), or 3) Full inertia-coordination model ("inertCoord") - each person's intercept, each person's own observed state variable lagged at the amount specified during the dataPrep step (again with linear trends removed), and each person's partner's state variable lagged at the amount specified (again with linear trends removed).
#'
#' @param basedata A dataframe that was produced with the "dataPrep" function.
#' @param whichModel Whether the model to be estimated is the inertia only model ("inert"), the coordination only model ("coord"), or the full inertia-coordination model ("inertCoord").
#' 
#' @return The function returns a list including: 1) the adjusted R^2 for the model for each dyad (called "R2") and 2) the parameter estimates for the model for each dyad (called "paramData", for use in either predicting, or being predicted by, the system variable).

#' @export

indivInertCoord <- function(basedata, whichModel)
{	
  if(whichModel != "inert" & whichModel != "coord" & whichModel != "inertCoord") {
  	stop("the model type must be either inert, coord or inertCoord")
	
	} else if (whichModel == "inert"){
	  model <- formula(obs_deTrend ~ -1 + dist0 + dist1 + dist0:obs_deTrend_Lag + dist1:obs_deTrend_Lag)
	  paramNames <- c("int0","int1","inert0","inert1","dyad") 
      
      } else if (whichModel == "coord"){
      	model <- formula(obs_deTrend ~ -1 + dist0 + dist1 + dist0:p_obs_deTrend_Lag + dist1:p_obs_deTrend_Lag)
      	paramNames <- c("int0","int1","coord0","coord1","dyad")
        
        } else {
          model <- formula(obs_deTrend ~ -1 + dist0 + dist1 + dist0:obs_deTrend_Lag + dist0:p_obs_deTrend_Lag + dist1:obs_deTrend_Lag + dist1:p_obs_deTrend_Lag)
          paramNames <- c("int0","int1","inert0","coord0","inert1","coord1","dyad")
  }
  
  newDiD <- unique(factor(basedata$dyad))
  R2 <- vector()
  param <- list()
	
  for (i in 1:length(newDiD)){  
    datai <- basedata[basedata$dyad == newDiD[i], ]
	m <- lm(model, na.action=na.exclude, data=datai)
	R2[[i]] <- summary(m)$adj.r.squared
	param[[i]] <- round(as.numeric(m$coefficients), 5)
	numParam <- length(m$coefficients)
	param[[i]][numParam + 1] <- unique(datai$dyad)
  }
  
  param <- as.data.frame(do.call(rbind, param))
  colnames(param) <- paramNames
  temp <- subset(basedata, select=c("id","dyad","sysVar","dist0"))
  temp2 <- unique(temp)
  paramData <- suppressMessages(plyr::join(param, temp2))
  
  results <- list(R2=R2, paramData=paramData)
}


#' Compares model fit for the inertia-only, coordination-only and full inertia-coordination model for each dyad's state trajectories using an R-square comparison. 
#' 
#' Fits inertia-only, coordination-only and full inertia-coordination models to each dyad's observed state variables and returns the adjusted R-squares, along with the differences between them, so positive values indicate better fit for the first model in the comparison. The 3 comparisons are inertia minus coordination, full model minus inertia, and full model minus coordination.
#'
#' @param basedata A dataframe that was produced with the "dataPrep" function.
#' 
#' @return The function returns a named list including: 1) the adjusted R^2 for the inertia model for each dyad (called "R2inert"), 2) the adjusted R^2 for the coordination model for each dyad (called "R2coord"), 3) the adjusted R^2 for the full inertia-coordination model for each dyad (called "R2inertCoord"), 4) the difference between the R-squares for each dyad for inertia minus coordination (called "R2dif_I_C"), 5) the difference for the full model minus inertia (called "R2dif_IC_I"), and 6) the difference for the full model minus coordination (called "R2dif_IC_C")

#' @export

indivInertCoordCompare <- function(basedata)
{ 
  newDiD <- unique(factor(basedata$dyad))
  R2inert <- vector()
  R2coord <- vector()
  R2inertCoord <- vector()
  R2dif_I_C <- vector()
  R2dif_IC_I <- vector()
  R2dif_IC_C <- vector()

  for (i in 1:length(newDiD)){
    datai <- basedata[basedata$dyad == newDiD[i], ]
	m1 <- formula(obs_deTrend ~ -1 + dist0 + dist1 + dist0:obs_deTrend_Lag + dist1:obs_deTrend_Lag)
	inert <- lm(m1, na.action=na.exclude, data=datai)
	m2 <- formula(obs_deTrend ~ -1 + dist0 + dist1 + dist0:p_obs_deTrend_Lag + dist1:p_obs_deTrend_Lag)
	coord <- lm(m2, na.action=na.exclude, data=datai)
    m3 <- formula(obs_deTrend ~ -1 + dist0 + dist1 + dist0:obs_deTrend_Lag + dist0:p_obs_deTrend_Lag + dist1:obs_deTrend_Lag + dist1:p_obs_deTrend_Lag)
    inertCoord <- lm(m3, na.action=na.exclude, data=datai)
			
	R2inert[[i]] <- summary(inert)$adj.r.squared
	R2coord[[i]] <- summary(coord)$adj.r.squared	
	R2inertCoord[[i]] <- summary(inertCoord)$adj.r.squared
			
	R2dif_I_C[[i]] <- R2inert[[i]] - R2coord[[i]]
	R2dif_IC_I[[i]] <- R2inertCoord[[i]] - R2inert[[i]]
	R2dif_IC_C[[i]] <- R2inertCoord[[i]] - R2coord[[i]]
  }			
  
  output <- list(R2inert=R2inert, R2coord=R2coord, R2inertCoord=R2inertCoord, R2dif_I_C=R2dif_I_C, R2dif_IC_I=R2dif_IC_I, R2dif_IC_C=R2dif_IC_C)
}


#' Produces plots of versions of the inertia-coordination model-predicted trajectories overlaid on raw data for each dyad.
#' 
#' The observed state variables (with linear trends removed) are predicted from one of the 3 versions of the inertia-coordination model (inertia only, "inert"; coordination only, "coord"; full inertia-coordination, "inertCoord") for each dyad individually. The predicted trajectories are plotted overlaid on the observed trajectories. 
#'
#' @param basedata A dataframe that was produced with the "dataPrep" function.
#' @param whichModel Whether the model to be estimated is the inertia only model ("inert"), the coordination only model ("coord"), or the full inertia-coordination model ("inertCoord").
#' @param dist0name A name for the level-0 of the distinguishing variable (e.g., "Women").
#' @param dist1name A name for the level-1 of the distinguishing variable (e.g., "Men").
#' @param obsName A name for the observed state variables being plotted (e.g., "Emotional Experience").
#' @param minMax An optional vector with desired minimum and maximum quantiles to be used for setting the y-axis range on the plots, e.g., minMax <- c(.1, .9) would set the y-axis limits to the 10th and 90th percentiles of the observed state variables. If not provided, the default is to use the minimum and maximum observed values of the state variables.
#' 
#' @return The function returns plots of the predicted values against the observed values for each dyad (called "plots"). The plots are also written to the working directory as a pdf file called "inertPlots.pdf", or "coordPlots.pdf" or "inertCoordPlots.pdf"

#' @import ggplot2
#' @export

indivInertCoordPlots <- function(basedata, whichModel, dist0name = NULL, dist1name = NULL, obsName = NULL, minMax=NULL)
{
  if(is.null(dist0name)){dist0name <- "dist0"}
  if(is.null(dist1name)){dist1name <- "dist1"}
  if(is.null(obsName)){obsName <- "observed"}

  if(whichModel != "inert" & whichModel != "coord" & whichModel != "inertCoord") {
  	stop("the model type must be either inert, coord or inertCoord")
	} else if (whichModel == "inert"){
	  model <- formula(obs_deTrend ~ -1 + dist0 + dist1 + dist0:obs_deTrend_Lag + dist1:obs_deTrend_Lag)
	  plotFileName <- "inertPlots.pdf"
      } else if (whichModel == "coord"){
      	model <- formula(obs_deTrend ~ -1 + dist0 + dist1 + dist0:p_obs_deTrend_Lag + dist1:p_obs_deTrend_Lag)
      	plotFileName <- "coordPlots.pdf"
        } else {
          model <- formula(obs_deTrend ~ -1 + dist0 + dist1 + dist0:obs_deTrend_Lag + dist0:p_obs_deTrend_Lag + dist1:obs_deTrend_Lag + dist1:p_obs_deTrend_Lag)
          plotFileName <- "inertCoordPlots.pdf"
        }

  if(is.null(minMax)){
  	min <- min(basedata$obs_deTrend, na.rm=T)
	max <- max(basedata$obs_deTrend, na.rm=T)
  } else {
  	min <- quantile(basedata$obs_deTrend, minMax[1], na.rm=T)
	max <- quantile(basedata$obs_deTrend, minMax[2],  na.rm=T)
  }

  newDiD <- unique(factor(basedata$dyad))
  plots <- list()
	
  for (i in 1:length(newDiD)){
	datai <- basedata[basedata$dyad == newDiD[i], ]
	m <- lm(model, na.action=na.exclude, data=datai)
	datai$obsPred <- predict(m)
	datai$role <- factor(datai$dist0, levels=c(0,1), labels=c(dist1name, dist0name)) 
	plotTitle <- as.character(unique(datai$dyad))
						
	plots[[i]] <- ggplot(datai, ggplot2::aes(x=time)) +
	geom_line(aes(y= obs_deTrend, color=role), linetype="dotted", size= .8, na.rm=T) +
	geom_line(aes(y=obsPred, color=role), size= .8, na.rm=T) + 
	scale_color_manual(name="Role", values=c("blue","red")) +
	ylab(obsName) +
	ylim(min, max) +
	annotate("text", x=-Inf, y=-Inf, hjust=0, vjust=0, label="Dots = Observed; Lines = Predicted", size=2) +
	labs(title= "Dyad ID:", subtitle= plotTitle) +
	theme(plot.title=element_text(size=11)) +
	theme(plot.subtitle=element_text(size=10))			
  }
	
  modelPlots <- gridExtra::marrangeGrob(grobs= plots, ncol=2, nrow=3)
  ggsave(plotFileName, modelPlots)
  results <- list(plots=plots)
}


#' Provides results for predicting the system variable from latent profiles of the inertia-coordination model parameters (obtained from latent profile analysis(LPA)). 
#' 
#' The system variable can be either dyadic (sysVarType = "dyadic"), where both partners have the same score (e.g., relationship length) or individual (sysVarType = "indiv"), where the partners can have different scores (e.g., age). For dyadic system variables, the only predictor is profile membership and the model is a regular regression model since all variables are at the level of the dyad. If the system variable is individual then the model is a random-intercept dyadic model and 3 models are estimated: 1) the main effect of profile membership, 2) main effects of profile membership and the distinguishing variable, and 3) the interaction of profile membership and the distinguishing variable. If the system variable is not normally distributed, any of the generalized linear models supported by glm (for dyadic system variables) or glmmPQL (for individual system variables) are available by specifying the "family" distribution.
#' 
#' @param basedata A dataframe containing the LPA profile memberships.
#' @param sysVarType Whether the system variable is "dyadic", which means both partners have the same score, or "indiv" which means the partners can have different scores
#' @param dist0name An optional name for the level-0 of the distinguishing variable (e.g., "Women"). Default is dist0.
#' @param dist1name An optional name for the level-1 of the distinguishing variable (e.g., "Men"). Default is dist1
#' @param sysVarName An optional name for the system variable being predicted (e.g., "Satisfaction"). Default is sysVar.
#' @param minMax An optional vector with desired minimum and maximum quantiles to be used for setting the y-axis range on the plots, e.g., minMax <- c(.1, .9) would set the y-axis limits to the 10th and 90th percentiles of the observed state variables. If not provided, the default is to use the minimum and maximum observed values of the state variables.
#' @param family An optional argument specifying the error distribution and link function to be used in the model. Any of the "family" options supported by glm (for dyadic system variables) or glmmPQL (for individual system variables) are available. Default is gaussian.
#' 
#' @return For normally distributed system variables, the function returns a list including the lm or lme objects containing the full results for each model (called "models"). Similarly, for non-normal system variables, the function returns a list of the glm or glmmPQL objects containing the full results for the models. By default, the function also displays histograms of the residuals and plots of the predicted values against observed values for each model, but these can be turned off by setting plots=F. 

#' @export
inertCoordSysVarOut <- function(basedata, sysVarType, dist0name=NULL, dist1name=NULL, sysVarName=NULL, minMax=NULL, family=NULL, printPlots=T)
{
  if(is.null(dist0name)){dist0name <- "dist0"}
  if(is.null(dist1name)){dist1name <- "dist1"}
  if(is.null(sysVarName)){sysVarName <- "sysVar"}
  if(is.null(family)){family <- "gaussian"}
	
  if(is.null(minMax)){
  	min <- min(basedata$sysVar, na.rm=T)
	max <- max(basedata$sysVar, na.rm=T)
  } else {
  	min <- quantile(basedata$sysVar, minMax[1], na.rm=T)
	max <- quantile(basedata$sysVar, minMax[2],  na.rm=T)
    }
 	
  if(sysVarType != "indiv" & sysVarType != "dyadic") {
	stop("the sysVarType must be either indiv or dyadic")
  }
	
  basedata <- basedata[complete.cases(basedata), ] 
  basedata$dist1 <- ifelse(basedata$dist0 == 1, 0, 1)
  basedata$dist <- factor(basedata$dist0, labels=c(dist1name, dist0name))
	
  if (sysVarType == "dyadic"){	
	basedata <- basedata[!duplicated(basedata$dyad), ]	
	
	if (family == "gaussian"){
	  profile <- lm(sysVar ~ profile, data= basedata)
	  profilePred <- fitted(profile)
	} else {
	  profile <- glm(sysVar ~ profile, data= basedata, family=family)
	  profilePred <- fitted(profile)
    }
      
    if(plots==T){
	  ylabName <- paste(sysVarName, "predicted", sep="_")
	  xlabName <- paste(sysVarName, "observed", sep="_")

	  hist(residuals(profile))
	  plot(profilePred ~ basedata$sysVar, xlim=c(min, max), ylim=c(min, max), ylab=ylabName, xlab=xlabName, main="Profile Model")
	 }
  }

  if (sysVarType == "indiv"){
	if (family == "gaussian"){
	
	profile <- nlme::lme(sysVar ~ profile, random= ~ 1 | dyad, data= basedata, na.action=na.omit, control=nlme::lmeControl(opt="optim"), method="ML")
	
	profilePlusDist <- nlme::lme(sysVar ~ profile + dist, random= ~ 1 | dyad, data= basedata, na.action=na.omit, control=nlme::lmeControl(opt="optim"), method="ML")

	profileByDist <- nlme::lme(sysVar ~ profile * dist, random= ~ 1 | dyad, data= basedata, na.action=na.omit, control=nlme::lmeControl(opt="optim"), method="ML")
    } else {
	
	  profile <- MASS::glmmPQL(sysVar ~ profile, random= ~ 1 | dyad, data= basedata, na.action=na.omit, control=nlme::lmeControl(opt="optim"), family=family)
	
	  profilePlusDist <- MASS::glmmPQL(sysVar ~ profile + dist, random= ~ 1 | dyad, data= basedata, na.action=na.omit, control=nlme::lmeControl(opt="optim"), family=family)

	  profileByDist <- MASS::glmmPQL(sysVar ~ profile * dist, random= ~ 1 | dyad, data= basedata, na.action=na.omit, control=nlme::lmeControl(opt="optim"), family=family)
	  }
	
	profilePred <- fitted(profile)
	profilePlusDistPred <- fitted(profilePlusDist)
	profileByDistPred <- fitted(profileByDist)
  	
    if(printPlots == T){
	  ylabName <- paste(sysVarName, "predicted", sep="_")
	  xlabName <- paste(sysVarName, "observed", sep="_")
	  
	  hist(residuals(profile))
	  plot(profilePred ~ basedata$sysVar, xlim=c(min, max), ylim=c(min, max), ylab=ylabName, xlab=xlabName, main="Profile Model")
	
	  hist(residuals(profilePlusDist))
	  plot(profilePlusDistPred ~ basedata$sysVar, xlim=c(min, max), ylim=c(min, max), ylab=ylabName, xlab=xlabName, main="Profile Plus Dist Model")

	  hist(residuals(profileByDist))
	  plot(profileByDistPred ~ basedata$sysVar, xlim=c(min, max), ylim=c(min, max), ylab=ylabName, xlab=xlabName, main="Profile By Dist Model")
	  
	  interaction.plot(basedata$profile, basedata$dist, basedata$sysVar, xlab="profile", ylab=sysVarName, trace.label="dist")
	}
  }

  if(sysVarType == "dyadic"){
	models <- list(profile=profile)
  }
	
  if(sysVarType == "indiv"){
	models <- list(profile=profile, profilePlusDist=profilePlusDist, profileByDist=profileByDist)
  }
	output <- list(models=models)
}


###############
#' Plots the bivariate state variables' model-predicted temporal trajectories for each latent profile of inertia-coordination parameters. 
#' 
#' Produces sets of prototypical example plots of the state variables' predicted temporal trajectories for each latent profile obtained based on the inertia-coordination parameters. The plots are produced by using the inertia-coordination parameters to predict temporal trajectories, with random noise added at each temporal step.  
#' 
#' @param origData A dataframe that was produced with the "dataPrep" function.
#' @param lpaData A dataframe containing the LPA profile memberships.
#' @param lpaParams A matrix containing the inertia-coordination parameter estimates associated with each of the latent profiles.
#' @param n_profiles The number of latent profiles.
#' @param time_length An optional value specifying how many time points to plot across. Default is 20.
#' @param dist0name An optional name for the level-0 of the distinguishing variable (e.g., "Women"). Default is dist0.
#' @param dist1name An optional name for the level-1 of the distinguishing variable (e.g., "Men"). Default is dist1
#' @param obsName An optional name for the state variables being plotted (e.g., "heart rate"). Default is obsName.
#' @param minMax An optional vector with desired minimum and maximum quantiles to be used for setting the y-axis range on the plots, e.g., minMax <- c(.1, .9) would set the y-axis limits to the 10th and 90th percentiles of the observed state variables. If not provided, the default is to use the minimum and maximum observed values of the state variables.
#' @param printPlots An optional argument specifying whether the plots should be automatically printed. Default is TRUE.
#' @param numPlots An optional value controlling how many random examples of each profile are produced. Default is 5.
#' 
#' @return The function returns the plots as a list. 

#' @import ggplot2
#' @export

inertCoordPredTraj <- function(origData, lpaData, lpaParams, n_profiles, time_length=NULL, dist0name=NULL, dist1name=NULL, obsName=NULL, minMax=NULL, printPlots=T, numPlots=NULL){

  if(is.null(time_length)){time_length <- 20}
  if(is.null(dist0name)){dist0name <- "dist0"}
  if(is.null(dist1name)){dist1name <- "dist1"}
  if(is.null(obsName)){obsName <- "observed"}
  
  if(is.null(minMax)){
  	min <- min(origData$obs_deTrend, na.rm=T)
	max <- max(origData$obs_deTrend, na.rm=T)
  } else {
  	min <- quantile(origData$obs_deTrend, minMax[1], na.rm=T)
	max <- quantile(origData$obs_deTrend, minMax[2],  na.rm=T)
  }

  noiseModel <- nlme::lme(obs_deTrend ~ -1 + dist0 + dist1 + dist0:obs_deTrend_Lag + dist0:p_obs_deTrend_Lag + dist1:obs_deTrend_Lag + dist1:p_obs_deTrend_Lag, random = ~ dist0 + dist1 | dyad, na.action=na.omit, data=origData, control=nlme::lmeControl(opt="optim"))
  
  noise <- noiseModel$sigma

  if(is.null(numPlots)) {numPlots <- 5}
  
  multiPlots <- list()
  plots <- list()
  label <- vector()

  
  for(i in 1:n_profiles){
  	for (k in 1:numPlots){
      start1 <- colMeans(subset(lpaData, profile==i & dist0==0, select=int1), na.rm=T)
      start0 <- colMeans(subset(lpaData, profile==i & dist0==0, select=int0), na.rm=T) 
      start <- c(start1, start0)
      paramsi <- lpaParams[ ,i]
      A <- matrix(paramsi, ncol=2, byrow=T)

      results1 <- list()
      results0 <- list()
      nextStep <- list()

      for (t in 1:time_length){
        if(t == 1){
  	      pred <- A %*% start
  	      dist <- c(1, 0)
  	      time <- t
  	      results1[[t]] <- list(pred=pred[1], dist=dist[1], time=time)
  	      results0[[t]] <- list(pred=pred[2], dist=dist[2], time=time)
   	      nextStep[[t]] <- pred + c(rnorm(n=2, mean=0, sd=noise))
        } else {
  	      pred <- A %*% nextStep[[t-1]]
  	      dist <- c(1, 0)
  	      time <- t
  	      results1[[t]] <- list(pred=pred[1], dist=dist[1], time=time)
  	      results0[[t]] <- list(pred=pred[2], dist=dist[2], time=time)
  	      nextStep[[t]] <- pred + c(rnorm(n=2, mean=0, sd=noise))
          }
        }
      
      final1 <- data.frame(do.call(rbind, results1))
      temp1 <- data.frame(matrix(unlist(final1), ncol=3, byrow=F))
      colnames(temp1) <- c("pred1", "dist1", "time")
      final0 <- data.frame(do.call(rbind, results0))
      temp0 <- data.frame(matrix(unlist(final0), ncol=3, byrow=F))
      colnames(temp0) <- c("pred0", "dist0", "time") 
      temp2 <- suppressMessages(plyr::join(temp1, temp0))
      temp3 <- reshape(temp2, idvar="dist", varying=list(c("pred1", "pred0"), c("dist1", "dist0")), direction="long")
      temp4 <- temp3[ , - 1]
      colnames(temp4) <- c("pred","dist","time")
      temp4$dist <- factor(temp4$dist, labels=c(dist0name, dist1name))
    
      plotData <- temp4
      profileName <- paste("Profile", i , sep="_")
      plots[[k]] <- ggplot(plotData, aes(x=time, y=pred, group=dist)) +
                  geom_line(aes(color=dist)) +
                  scale_color_manual(values=c("black","gray47")) +
                  ylab(obsName) +
                  ylim(min, max) +
	              labs(title= profileName, subtitle= "Predicted_Trajectory") +
	              theme(plot.title=element_text(size=11)) 
    }
    multiPlots[[i]] <- plots
  }
  if(printPlots == T) {print(multiPlots)}
  return(multiPlots)
}










