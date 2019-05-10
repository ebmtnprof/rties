

######## This file includes all the functions needed for an inertia-coordination analysis

################### autoCorPlots

#' Produces auto-correlation plots of the observed state variable for lags of -+ 20 time steps for each dyad.
#' 
#' @param basedata A user provided dataframe.
#' @param dyadId The name of the column in the dataframe that has the dyad-level identifier.
#' @param personId The name of the column in the dataframe that has the person-level identifier.
#' @param obs_name The name of the column in the dataframe that has the time-varying observable (e.g., the variable for which dynamics will be assessed).
#' @param time_name The name of the column in the dataframe that indicates sequential temporal observations.

#' @export

autoCorPlots <- function(basedata, dyadId, personId, obs_name, time_name){invisible(acp(basedata, personId, dyadId, obs_name, time_name))}

acp <- function(basedata, personId, dyadId, obs_name, time_name)
{
  basedata <- subset(basedata, select=c(personId, dyadId, obs_name, time_name))
  names(basedata)[1] <- "id"
  names(basedata)[2] <- "dyad"
  names(basedata)[3] <- "obs"
  names(basedata)[4] <- "time"

  basedata <- basedata[complete.cases(basedata), ]
  basedata <- lineCenterById(basedata)
  basedata <- actorPartnerDataTime(basedata, "dyad", "id")
  
  newDiD <- unique(factor(basedata$dyad))
  acf <- list()
  plotTitle <- list()
	
  for (i in 1:length(newDiD)){
	datai <- basedata[basedata$dyad == newDiD[i], ]
	datai_ts <- zoo::zoo(datai[,"obs_deTrend"])
	d <- acf(datai_ts, plot=F, na.action=na.exclude, lag.max=20)
    acf[[i]] <- d
    plotTitle[[i]] <- as.character(unique(datai$dyad))
  }  
  
  par(mfrow=c(2,2))
  for(j in 1:length(newDiD)){
  plot(acf[[j]], main=paste("AutoCorr_Dyad", plotTitle[j], sep="_"))
  }
  par(mfrow=c(1,1))
  return(acf)
}

###################### crossCorPlots

#' Produces cross-correlation plots of the observed state variable for lags of -+ 20 time steps for each dyad.
#' 
#' @param basedata A user provided dataframe.
#' @param dyadId The name of the column in the dataframe that has the dyad-level identifier.
#' @param personId The name of the column in the dataframe that has the person-level identifier.
#' @param obs_name The name of the column in the dataframe that has the time-varying observable (e.g., the variable for which dynamics will be assessed).
#' @param time_name The name of the column in the dataframe that indicates sequential temporal observations.

#' @export

crossCorPlots <- function(basedata, personId, dyadId, obs_name, time_name){invisible(ccp(basedata, personId, dyadId, obs_name, time_name))}

ccp <- function(basedata, personId, dyadId, obs_name, time_name)
{
  basedata <- subset(basedata, select=c(personId, dyadId, obs_name, time_name))  
  names(basedata)[1] <- "id"
  names(basedata)[2] <- "dyad"
  names(basedata)[3] <- "obs"
  names(basedata)[4] <- "time"

  basedata <- basedata[complete.cases(basedata), ]
  basedata <- lineCenterById(basedata)
  basedata <- actorPartnerDataTime(basedata, "dyad", "id")
  
  newDiD <- unique(factor(basedata$dyad))
  ccf <- list()
  plotTitle <- list()
	
  for (i in 1:length(newDiD)){
	datai <- basedata[basedata$dyad == newDiD[i], ]
	datai_ts1 <- zoo::zoo(datai[ ,"obs_deTrend"])
	datai_ts2 <- zoo::zoo(datai[ ,"p_obs_deTrend"])
	d <- ccf(datai_ts1, datai_ts2, type="correlation", plot=F, na.action=na.exclude, lag.max=20)
    ccf[[i]] <- d
    plotTitle[[i]] <- as.character(unique(datai$dyad))
  }  
  
  par(mfrow=c(2,2))
  for(j in 1:length(newDiD)){
  plot(ccf[[j]], main=paste("CrossCorr_Dyad", plotTitle[j], sep="_"))
  }
  par(mfrow=c(1,1))
  return(acf)
}


######################## indivInertCoord

#' Estimates versions of the inertia-coordination model for each dyad.
#' 
#' The user specifies which of 3 models are to be estimated. Each model predicts the observed state variables (with linear trends removed) from either: 1) Inertia only ("inert")- each person's intercept and each person's own observed state variable lagged at the amount specified during the dataPrep step (again with linear trends removed), 2) Coordination only ("coord")- each person's intercept and each person's partner's state variable lagged at the amount specified (again with linear trends removed), or 3) Full inertia-coordination model ("inertCoord") - each person's intercept, each person's own observed state variable lagged at the amount specified during the dataPrep step (again with linear trends removed), and each person's partner's state variable lagged at the amount specified (again with linear trends removed).
#'
#' @param prepData A dataframe that was produced with the "dataPrep" function.
#' @param whichModel Whether the model to be estimated is the inertia only model ("inert"), the coordination only model ("coord"), or the full inertia-coordination model ("inertCoord").
#' 
#' @return The function returns a dataframe containing the parameter estimates, called "params", for use in the latent profile analysis.

#' @export

indivInertCoord <- function(prepData, whichModel)
{	
  basedata <- prepData
  
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
  temp <- subset(basedata, select=c("id","dyad","dist0"))
  temp2 <- unique(temp)
  data <- suppressMessages(plyr::join(param, temp2))
  params <- data[data$dist0 == 1, ]
  
  results <- list(params=params)
}


####################### indivInertCoordCompare

#' Compares model fit for the inertia-only, coordination-only and full inertia-coordination model for each dyad's state trajectories using an R-square comparison. 
#' 
#' Fits inertia-only, coordination-only and full inertia-coordination models to each dyad's observed state variables and returns the adjusted R-squares, along with the differences between them, so positive values indicate better fit for the first model in the comparison. The 3 comparisons are inertia minus coordination, full model minus inertia, and full model minus coordination.
#'
#' @param prepData A dataframe that was produced with the "dataPrep" function.
#' 
#' @return The function returns a named list including: 1) the adjusted R^2 for the inertia model for each dyad (called "R2inert"), 2) the adjusted R^2 for the coordination model for each dyad (called "R2coord"), 3) the adjusted R^2 for the full inertia-coordination model for each dyad (called "R2inertCoord"), 4) the difference between the R-squares for each dyad for inertia minus coordination (called "R2dif_I_C"), 5) the difference for the full model minus inertia (called "R2dif_IC_I"), and 6) the difference for the full model minus coordination (called "R2dif_IC_C")

#' @export

indivInertCoordCompare <- function(prepData)
{ 
  basedata <- prepData
  
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


####################### indivInertCoordPlots

#' Produces plots of the inertia-coordination model-predicted trajectories overlaid on raw data for each dyad.
#' 
#' The observed state variables (with linear trends removed) are predicted from one of the 3 versions of the inertia-coordination model (inertia only, "inert"; coordination only, "coord"; full inertia-coordination, "inertCoord") for each dyad individually. The predicted trajectories are plotted overlaid on the observed trajectories. 
#'
#' @param prepData A dataframe that was produced with the "dataPrep" function.
#' @param whichModel Whether the model to be estimated is the inertia only model ("inert"), the coordination only model ("coord"), or the full inertia-coordination model ("inertCoord").
#' @param dist0name An optional name for the level-0 of the distinguishing variable to appear on plots (e.g., "Women").
#' @param dist1name An optional name for the level-1 of the distinguishing variable to appear on plots (e.g., "Men").
#' @param plot_obs_name An optional name for the observed state variable to appear on plots (e.g., "Emotional Experience").
#' @param minMax An optional vector with desired minimum and maximum quantiles to be used for setting the y-axis range on the plots, e.g., minMax <- c(.1, .9) would set the y-axis limits to the 10th and 90th percentiles of the observed state variables. If not provided, the default is to use the minimum and maximum observed values of the state variables.
#' 
#' @return The function returns plots of the predicted values against the observed values for each dyad (called "plots"). The plots are also written to the working directory as a pdf file called "inertPlots.pdf", or "coordPlots.pdf" or "inertCoordPlots.pdf"

#' @import ggplot2
#' @export

indivInertCoordPlots <- function(prepData, whichModel, dist0name = NULL, dist1name = NULL, plot_obs_name = NULL, minMax=NULL)
{
  basedata <- prepData
  
  if(is.null(dist0name)){dist0name <- "dist0"}
  if(is.null(dist1name)){dist1name <- "dist1"}
  if(is.null(plot_obs_name)){plot_obs_name <- "observed"}

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
	ylab(plot_obs_name) +
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

#################### inertCoordResids

#' Produces histograms of the residuals from the inertia-coordination model for each dyad.
#' 
#' @param prepData A dataframe that was produced with the "dataPrep" function.
#' @param whichModel Whether the model to be estimated is the inertia only model ("inert"), the coordination only model ("coord"), or the full inertia-coordination model ("inertCoord").
#' 
#' @return The function returns histograms of the residuals from the model for each dyad (called "plots"). The plots are also written to the working directory as a pdf file called "inertResid.pdf", or "coordResid.pdf" or "inertCoordResid.pdf"

#' @import ggplot2
#' @export

inertCoordResids <- function(prepData, whichModel)
{
  basedata <- prepData
  
  if(whichModel != "inert" & whichModel != "coord" & whichModel != "inertCoord") {
  	stop("the model type must be either inert, coord or inertCoord")
	} else if (whichModel == "inert"){
	  model <- formula(obs_deTrend ~ -1 + dist0 + dist1 + dist0:obs_deTrend_Lag + dist1:obs_deTrend_Lag)
	  plotFileName <- "inertResid.pdf"
      } else if (whichModel == "coord"){
      	model <- formula(obs_deTrend ~ -1 + dist0 + dist1 + dist0:p_obs_deTrend_Lag + dist1:p_obs_deTrend_Lag)
      	plotFileName <- "coordResid.pdf"
        } else {
          model <- formula(obs_deTrend ~ -1 + dist0 + dist1 + dist0:obs_deTrend_Lag + dist0:p_obs_deTrend_Lag + dist1:obs_deTrend_Lag + dist1:p_obs_deTrend_Lag)
          plotFileName <- "inertCoordResid.pdf"
        }

  newDiD <- unique(factor(basedata$dyad))
  plots <- list()
  resid <- list()
	
  for (i in 1:length(newDiD)){
	datai <- basedata[basedata$dyad == newDiD[i], ]
	m <- lm(model, na.action=na.exclude, data=datai) 
	plotTitle <- as.character(unique(datai$dyad))
	resid[[i]] <- m$residuals
	plotResid <- data.frame(resid[[i]])
	colnames(plotResid) <- "Residuals"
						
	plots[[i]] <- ggplot(plotResid, aes(x=Residuals)) +
	geom_histogram(color="black", fill="grey") +
	labs(title= "Dyad ID:", subtitle= plotTitle) +
	theme(plot.title=element_text(size=11)) +
	theme(plot.subtitle=element_text(size=10))		
  }
	
  modelPlots <- suppressMessages(gridExtra::marrangeGrob(grobs= plots, ncol=2, nrow=3))
  ggsave(plotFileName, modelPlots)
  results <- list(plots=plots)
}


######################## inertCoordPlotTrajInternal 

#' Plots the bivariate state variables' model-predicted temporal trajectories for each latent profile of inertia-coordination parameters. Not exported - used by inspectProfiles
#' 
#' Produces sets of prototypical example plots of the state variables' predicted temporal trajectories for each latent profile obtained based on the inertia-coordination parameters. The plots are produced by using the inertia-coordination parameters to predict temporal trajectories, with random noise added at each temporal step.  
#' 
#' @param prepData A dataframe that was produced with the "dataPrep" function.
#' @param paramEst A dataframe created by indivInertCoord containing the inertia-coordination parameter estimates for each dyad.
#' @param n_profiles The number of latent profiles.
#' @param dist0name An optional name for the level-0 of the distinguishing variable (e.g., "Women"). Default is dist0.
#' @param dist1name An optional name for the level-1 of the distinguishing variable (e.g., "Men"). Default is dist1
#' @param minMax An optional vector with desired minimum and maximum quantiles to be used for setting the y-axis range on the plots, e.g., minMax <- c(.1, .9) would set the y-axis limits to the 10th and 90th percentiles of the observed state variables. If not provided, the default is to use the minimum and maximum observed values of the state variables.
#' 
#' @return The function prints the plots. 

#' @import ggplot2
#' @import dplyr

inertCoordPlotTrajInternal <- function(prepData, paramEst, n_profiles, time_lag, dist0name=NULL, dist1name=NULL, minMax=NULL)
{ 
  time_length <- 20
  if(is.null(dist0name)){dist0name <- "dist0"}
  if(is.null(dist1name)){dist1name <- "dist1"}
  
  if(is.null(minMax)){
  	min <- min(prepData$obs_deTrend, na.rm=T)
	max <- max(prepData$obs_deTrend, na.rm=T)
  } else {
  	min <- quantile(prepData$obs_deTrend, minMax[1], na.rm=T)
	max <- quantile(prepData$obs_deTrend, minMax[2],  na.rm=T)
  }

  temp1 <- paramEst %>%
        select(inert1, coord1, coord0, inert0) %>%
        tidyLPA::estimate_profiles(n_profiles)
        profileParams <- as.data.frame(tidyLPA::get_estimates(temp1))
  
  noiseModel <- nlme::lme(obs_deTrend ~ -1 + dist0 + dist1 + dist0:obs_deTrend_Lag + dist0:p_obs_deTrend_Lag + dist1:obs_deTrend_Lag + dist1:p_obs_deTrend_Lag, random = ~ dist0 + dist1 | dyad, na.action=na.exclude, data=prepData, control=nlme::lmeControl(opt="optim"))
  
  noise <- noiseModel$sigma

  numPlots <- 3
  
  stateCutOff <- time_lag + 5
  
  multiPlots <- list()
  plots <- list()
  label <- vector()

  for(i in 1:n_profiles){
  	for (k in 1:numPlots){
      statedata0 <- prepData[prepData$dist0 == 1 & prepData$time <= stateCutOff,] 
	  start0 <- median(statedata0$obs_deTrend, na.rm=T)
  	  statedata1 <- prepData[prepData$dist0 == 0 & prepData$time <= stateCutOff,] 
	  start1 <- median(statedata1$obs_deTrend, na.rm=T)
  	
      start <- c(start1, start0)  
      
      tempParams <- subset(profileParams, Category == "Means" & Class==i, select=c(Parameter, Estimate))
	  temp1 <- tempParams[ ,2]
      names <- tempParams[ ,1]
      names(temp1) <- names
      paramsi <- temp1

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
                  ylab("observed") +
                  ylim(min, max) +
	              labs(title= profileName, subtitle= "Predicted_Trajectory") +
	              theme(plot.title=element_text(size=11)) 
    }
    multiPlots[[i]] <- plots
  }
print(multiPlots)
}


######################## inertCoordPlotTraj 

#' Plots the bivariate state variables' model-predicted temporal trajectories for each latent profile of inertia-coordination parameters. 
#' 
#' Produces sets of prototypical example plots of the state variables' predicted temporal trajectories for each latent profile obtained based on the inertia-coordination parameters. The plots are produced by using the inertia-coordination parameters to predict temporal trajectories, with random noise added at each temporal step.  
#' 
#' @param prepData A dataframe that was produced with the "dataPrep" function.
#' @param paramEst A dataframe created by indivInertCoord containing the inertia-coordination parameter estimates for each dyad.
#' @param n_profiles The number of latent profiles.
#' @param time_length An optional value specifying how many time points to plot across. Default is 20.
#' @param dist0name An optional name for the level-0 of the distinguishing variable (e.g., "Women"). Default is dist0.
#' @param dist1name An optional name for the level-1 of the distinguishing variable (e.g., "Men"). Default is dist1
#' @param minMax An optional vector with desired minimum and maximum quantiles to be used for setting the y-axis range on the plots, e.g., minMax <- c(.1, .9) would set the y-axis limits to the 10th and 90th percentiles of the observed state variables. If not provided, the default is to use the minimum and maximum observed values of the state variables.
#' @param numPlots An optional value controlling how many random examples of each profile are produced. Default is 5.
#' 
#' @return The function prints the plots and returns them as a list. 

#' @import ggplot2
#' @import dplyr

#' @export

inertCoordPlotTraj <- function(prepData, paramEst, n_profiles, time_lag, time_length=NULL, dist0name=NULL, dist1name=NULL, minMax=NULL, numPlots=NULL)
{ 
  if(is.null(time_length)){time_length <- 20}
  if(is.null(dist0name)){dist0name <- "dist0"}
  if(is.null(dist1name)){dist1name <- "dist1"}
  
  if(is.null(minMax)){
  	min <- min(prepData$obs_deTrend, na.rm=T)
	max <- max(prepData$obs_deTrend, na.rm=T)
  } else {
  	min <- quantile(prepData$obs_deTrend, minMax[1], na.rm=T)
	max <- quantile(prepData$obs_deTrend, minMax[2],  na.rm=T)
  }

  temp1 <- paramEst %>%
        select(inert1, coord1, coord0, inert0) %>%
        tidyLPA::estimate_profiles(n_profiles)
        profileParams <- as.data.frame(tidyLPA::get_estimates(temp1))
  
  noiseModel <- nlme::lme(obs_deTrend ~ -1 + dist0 + dist1 + dist0:obs_deTrend_Lag + dist0:p_obs_deTrend_Lag + dist1:obs_deTrend_Lag + dist1:p_obs_deTrend_Lag, random = ~ dist0 + dist1 | dyad, na.action=na.omit, data=prepData, control=nlme::lmeControl(opt="optim"))
  
  noise <- noiseModel$sigma
  
  stateCutOff <- time_lag + 5

  if(is.null(numPlots)){numPlots <- 5}
  
  multiPlots <- list()
  plots <- list()
  label <- vector()

  for(i in 1:n_profiles){
  	for (k in 1:numPlots){
      statedata0 <- prepData[prepData$dist0 == 1 & prepData$time <= stateCutOff,] 
	  start0 <- median(statedata0$obs_deTrend, na.rm=T)
  	  statedata1 <- prepData[prepData$dist0 == 0 & prepData$time <= stateCutOff,] 
	  start1 <- median(statedata1$obs_deTrend, na.rm=T)
  	
      start <- c(start1, start0)  
      
      tempParams <- subset(profileParams, Category == "Means" & Class==i, select=c(Parameter, Estimate))
	  temp1 <- tempParams[ ,2]
      names <- tempParams[ ,1]
      names(temp1) <- names
      paramsi <- temp1

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
                  ylab("observed") +
                  ylim(min, max) +
	              labs(title= profileName, subtitle= "Predicted_Trajectory") +
	              theme(plot.title=element_text(size=11)) 
    }
    multiPlots[[i]] <- plots
  }
print(multiPlots)
return(multiPlots)
}









