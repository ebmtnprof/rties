

######## This file includes all the functions needed for an inertia-coordination analysis

#' Estimates an "inertia only" model for each dyad.
#' 
#' The observed state variables (with linear trends removed) are predicted from each person's intercept and each person's own observed state variable lagged at the amount specified during the dataPrep step (again with linear trends removed)
#'
#' @param basedata A dataframe that was produced with the "dataPrep" function.
#' @param dist0name A name for the level-0 of the distinguishing variable (e.g., "Women").
#' @param dist1name A name for the level-1 of the distinguishing variable (e.g., "Men").
#' @param obsName A name for the observed state variables being plotted (e.g., "Emotional Experience").
#' 
#' @return The function returns a list including: 1) the adjusted R^2 for the model for each dyad (called "R2"), 2) the parameter estimates for the model for each dyad (called "paramData", for use in either predicting, or being predicted by, the system variable), and 3) plots of the predicted values against the observed values for each dyad (called "plots"). The plots are also written to the working directory as a pdf file called "inertPlots.pdf"

#' @import ggplot2
#' @export
indivInert <- function(basedata, dist0name, dist1name, obsName)
{
	newDiD <- unique(factor(basedata$dyad))
	min <- min(basedata$obs_deTrend, na.rm=T)
	max <- max(basedata$obs_deTrend, na.rm=T)
	
	r2 <- vector()
	param <- list()
	plots <- list()
	
		for (i in 1:length(newDiD))
		{
			datai <- basedata[basedata$dyad == newDiD[i], ]
			m <- lm(obs_deTrend ~ -1 + dist0 + dist1 + dist0:obs_deTrend_Lag + dist1:obs_deTrend_Lag, na.action=na.exclude, data=datai)
			r2[[i]] <- summary(m)$adj.r.squared
			param[[i]] <- round(as.numeric(m$coefficients), 5)
			numParam <- length(m$coefficients)
			param[[i]][numParam + 1] <- unique(datai$dyad)
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
	
		param <- as.data.frame(do.call(rbind, param))
		colnames(param) <- c("int0","int1","inert0","inert1","dyad")
		temp <- subset(basedata, select=c("id","dyad","sysVar","dist0"))
		temp2 <- unique(temp)
		paramData <- plyr::join(param, temp2)
		
		inertPlots <- gridExtra::marrangeGrob(grobs= plots, ncol=2, nrow=3)
		ggsave('inertPlots.pdf', inertPlots)

	results <- list(R2=r2, paramData=paramData, plots=plots)
}


#' Estimates a "coordination only" model for each dyad.
#' 
#' The observed state variables (with linear trends removed) are predicted from each person's intercept and each person's partner's concurrent state variable (again with linear trends removed).
#' 
#' @param basedata A dataframe that was produced with the "dataPrep" function.
#' @param dist0name A name for the level-0 of the distinguishing variable (e.g., "Women").
#' @param dist1name A name for the level-1 of the distinguishing variable (e.g., "Men").
#' @param obsName A name for the observed state variables being plotted (e.g., "Emotional Experience").
#' 
#' @return The function returns a list including: 1) the adjusted R^2 for the model for each dyad (called "R2"), 2) the parameter estimates for the model for each dyad (called "paramData", for use in either predicting, or being predicted by, the system variable), and 3) plots of the predicted values against the observed values for each dyad (called "plots"). The plots are also written to the working directory as a pdf file called "coordPlots.pdf"

#' @import ggplot2
#' @export
indivCoord <- function(basedata, dist0name, dist1name, obsName)
{
	newDiD <- unique(factor(basedata$dyad))
	min <- min(basedata$obs_deTrend, na.rm=T)
	max <- max(basedata$obs_deTrend, na.rm=T)
	
	r2 <- vector()
	param <- list()
	plots <- list()
	
		for (i in 1:length(newDiD))
		{
			datai <- basedata[basedata$dyad == newDiD[i], ]
			m <- lm(obs_deTrend ~ -1 + dist0 + dist1 + dist0:p_obs_deTrend + dist1:p_obs_deTrend, na.action=na.exclude, data=datai)
			r2[[i]] <- summary(m)$adj.r.squared
			param[[i]] <- round(as.numeric(m$coefficients), 5)
			numParam <- length(m$coefficients)
			param[[i]][numParam + 1] <- unique(datai$dyad)
			datai$obsPred <- predict(m)
			datai$role <- factor(datai$dist0, levels=c(0,1), labels=c(dist1name, dist0name)) 
			plotTitle <- as.character(unique(datai$dyad))
						
			plots[[i]] <- ggplot(datai, aes(x=time)) +
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
	
		param <- as.data.frame(do.call(rbind, param))
		colnames(param) <- c("int0","int1","coord0","coord1","dyad")
		temp <- subset(basedata, select=c("id","dyad","sysVar","dist0"))
		temp2 <- unique(temp)
		paramData <- plyr::join(param, temp2)
		
		coordPlots <- gridExtra::marrangeGrob(grobs= plots, ncol=2, nrow=3)
		ggsave('coordPlots.pdf', coordPlots)

	results <- list(R2=r2, paramData=paramData, plots=plots)
}


#' Estimates the full "inertia-coordination" model for each dyad.
#' 
#' The observed state variables (with linear trends removed) are predicted from each person's intercept, each person's own observed state variable lagged at the amount specified during the dataPrep step (again with linear trends removed), and each person's partner's concurrent state variable (again with linear trends removed).
#' 
#' @param basedata A dataframe that was produced with the "dataPrep" function.
#' @param dist0name A name for the level-0 of the distinguishing variable (e.g., "Women").
#' @param dist1name A name for the level-1 of the distinguishing variable (e.g., "Men").
#' @param obsName A name for the observed state variables being plotted (e.g., "Emotional Experience").
#' 
#' @return The function returns a list including: 1) the adjusted R^2 for the model for each dyad (called "R2"), 2) the parameter estimates for the model for each dyad (called "paramData", for use in either predicting, or being predicted by, the system variable), and 3) plots of the predicted values against the observed values for each dyad (called "plots"). The plots are also written to the working directory as a pdf file called "inertCoordPlots.pdf"

#' @import ggplot2
#' @export
indivInertCoord <- function(basedata, dist0name, dist1name, obsName)
{
	newDiD <- unique(factor(basedata$dyad))
	min <- min(basedata$obs_deTrend, na.rm=T)
	max <- max(basedata$obs_deTrend, na.rm=T)
	
	r2 <- vector()
	param <- list()
	plots <- list()
	
		for (i in 1:length(newDiD))
		{
			datai <- basedata[basedata$dyad == newDiD[i], ]
			m <- lm(obs_deTrend ~ -1 + dist0 + dist1 + dist0:obs_deTrend_Lag + dist0:p_obs_deTrend + dist1:obs_deTrend_Lag + dist1:p_obs_deTrend, na.action=na.exclude, data=datai)
			r2[[i]] <- summary(m)$adj.r.squared
			param[[i]] <- round(as.numeric(m$coefficients), 5)
			numParam <- length(m$coefficients)
			param[[i]][numParam + 1] <- unique(datai$dyad)
			datai$obsPred <- predict(m)
			datai$role <- factor(datai$dist0, levels=c(0,1), labels=c(dist1name, dist0name)) 
			plotTitle <- as.character(unique(datai$dyad))
						
			plots[[i]] <- ggplot(datai, aes(x=time)) +
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
	
		param <- as.data.frame(do.call(rbind, param))
		colnames(param) <- c("int0","int1","inert0","coord0","inert1","coord1","dyad")
		temp <- subset(basedata, select=c("id","dyad","sysVar","dist0"))
		temp2 <- unique(temp)
		paramData <- plyr::join(param, temp2)
		
		inertCoordPlots <- gridExtra::marrangeGrob(grobs= plots, ncol=2, nrow=3)
		ggsave('inertCoordPlots.pdf', inertCoordPlots)

	results <- list(R2=r2, paramData=paramData, plots=plots)
}

#' Compares model fit for the inertia-only, coordination-only and full inertia-coordination model for each dyad's state trajectories using an R-square comparison. 
#' 
#' Fits inertia-only, coordination-only and full inertia-coordination models to each dyad's observed state variables and returns the adjusted R-squares, along with the differences between them, so positive values indicate better fit for the first model in the comparison. The 3 comparisons are inertia minus coordination, full model minus inertia, and full model minus coordination.
#'
#' @param basedata A dataframe that was produced with the "dataPrep" function.
#' 
#' @return The function returns a named list including: 1) the adjusted R^2 for the inertia model for each dyad (called "R2inert"), 2) the adjusted R^2 for the coordination model for each dyad (called "R2coord"), 3) the adjusted R^2 for the full inertia-coordination model for each dyad (called "R2inertCoord"), 4) the difference between the R-squares for each dyad for inertia minus coordination (called "R2dif_I_C"), 5) the difference for the full model minus inertia (called "R2dif_IC_I"), and 6) the difference for the full model minus coordination (called "R2dif_IC_C")
#' @export

inertCoordIndivCompare <- function(basedata)
{
	newDiD <- unique(factor(basedata$dyad))
	
	R2inert <- vector()
	R2coord <- vector()
	R2inertCoord <- vector()
	R2dif_I_C <- vector()
	R2dif_IC_I <- vector()
	R2dif_IC_C <- vector()

		for (i in 1:length(newDiD))
		{
			datai <- basedata[basedata$dyad == newDiD[i], ]
			
			inert <- lm(obs_deTrend ~ -1 + dist0 + dist1 + dist0:obs_deTrend_Lag + dist1:obs_deTrend_Lag, na.action=na.exclude, data=datai)
	
			coord <- lm(obs_deTrend ~ -1 + dist0 + dist1 + dist0:p_obs_deTrend + dist1:p_obs_deTrend, na.action=na.exclude, data=datai)

			inertCoord <- lm(obs_deTrend ~ -1 + dist0 + dist1 + dist0:obs_deTrend_Lag + dist0:p_obs_deTrend + dist1:obs_deTrend_Lag + dist1:p_obs_deTrend, na.action=na.exclude, data=datai)
			
			R2inert[[i]] <- summary(inert)$adj.r.squared
			R2coord[[i]] <- summary(coord)$adj.r.squared	
			R2inertCoord[[i]] <- summary(inertCoord)$adj.r.squared
			
			R2dif_I_C[[i]] <- R2inert[[i]] - R2coord[[i]]
			R2dif_IC_I[[i]] <- R2inertCoord[[i]] - R2inert[[i]]
			R2dif_IC_C[[i]] <- R2inertCoord[[i]] - R2coord[[i]]
		}			
		output <- list(R2inert=R2inert, R2coord=R2coord, R2inertCoord=R2inertCoord, R2dif_I_C=R2dif_I_C, R2dif_IC_I=R2dif_IC_I, R2dif_IC_C=R2dif_IC_C)
		}

#' Compares variants of the inertia-coordination model for predicting the system variable from the dynamic parameters.
#' 
#' The dynamic parameters used in these models come from a set including both people's inertia (inert0 and inert1) and coordination (coord0 and coord1) estimates. The 4 models compared are a baseline intercept-only, inertia-only (inert0 + inert1), coordination-only (coord0 + coord1) and full inertia-coordination (inert0 + inert1 + coord0 + coord1) models. The system variable can be either dyadic (sysVarType = "dyad"), where both partners have the same score (e.g., relationship length) or individual (sysVarType = "indiv"), where the partners can have different scores (e.g., age). If it is individual then both actor and partner effects of the dynamic parameters are included. 
#' 
#' @param basedata A dataframe containing the inertia-coordination parameter estimates produced by the "indivInertCoord" function.
#' @param sysVarType Whether the system variable is "dyad", which means both partners have the same socre, or "indiv" which means the partners can have different scores
#' @param dist0name A name for the level-0 of the distinguishing variable (e.g., "Women").
#' @param dist1name A name for the level-1 of the distinguishing variable (e.g., "Men").
#' @param sysVarName A name for the system variable being predicted (e.g., "Satisfaction").
#' 
#' @return The function returns a list including: 1) the lm or lme objects containing the full results for each model (called "models"), and 2) adjusted R^2 information for each model (called "R2"). The function also displays histograms of the residuals and plots of the predicted values against observed values for each model. 

#' @export
inertCoordSysVarOut <- function(basedata, sysVarType, dist0name, dist1name, sysVarName)
{
	basedata <- basedata[complete.cases(basedata), ] # remove NAs in order to compare nested models	
	basedata$dist1 <- ifelse(basedata$dist0 == 1, 0, 1)

	# Names for indiv model parameters
	aInert0name <- paste("actorInert",dist0name, sep="_")
	pInert0name <- paste("partnerInert",dist0name, sep="_")
	aInert1name <- paste("actorInert",dist1name, sep="_")
	pInert1name <- paste("partnerInert",dist1name, sep="_")
	aCoord0name <- paste("actorCoord",dist0name, sep="_")
	aCoord1name <- paste("actorCoord",dist1name, sep="_")
	pCoord0name <- paste("partnerCoord",dist0name, sep="_")
	pCoord1name <- paste("partnerCoord",dist1name, sep="_")
	
	# Names for dyad model parameters
	inert0name <- paste("inert",dist0name, sep="_")
	inert1name <- paste("inert",dist1name, sep="_")
	coord0name <- paste("coord",dist0name, sep="_")
	coord1name <- paste("coord",dist1name, sep="_")

 	if(sysVarType != "indiv" & sysVarType != "dyad") 
    {
	stop("the sysVarType must be either indiv or dyad")
	}
	
 	else if (sysVarType == "dyad")
 	{	
	#basedata <- basedata[!duplicated(basedata$dyad), ]
	
	# Base dyadic sysVar
	base <- lm(sysVar ~ 1, data= basedata)
	names(base$coefficients) <- c("intercept")

	# Inert dyadic sysVar 
	inert <- lm(sysVar ~ inert0 + inert1, data= basedata)
	names(inert$coefficients) <- c("intercept", inert0name, inert1name)

	# Coord dyadic sysVar 
	coord <- lm(sysVar ~ coord0 + coord1, data= basedata)
	names(coord$coefficients) <- c("intercept", coord0name, coord1name)

	# InertCoord dyadic sysVar 
	inertCoord <- lm(sysVar ~ inert0 + inert1 + coord0 + coord1, data= basedata)
	names(inertCoord$coefficients) <- c("intercept", inert0name, inert1name, coord0name, coord1name)

	basePred <- predict(base)## all predictions will be the mean
	inertPred <- predict(inert)
	coordPred <- predict(coord)
	inertCoordPred <- predict(inertCoord)

	baseR2 <- summary(base)$adj.r.squared # will be zero
	inertR2 <- summary(inert)$adj.r.squared
	coordR2 <- summary(coord)$adj.r.squared
	inertCoordR2 <- summary(inertCoord)$adj.r.squared
	}

	else if (sysVarType == "indiv")
	{
	# Base indiv sysVar
	base <- nlme::lme(sysVar ~ dist0, random= ~ 1 | dyad, data= basedata, na.action=na.omit, method="ML", control=nlme::lmeControl(opt="optim"))
	names(base$coefficients$fixed) <- c("intercept", dist0name)
	
	# Inert indiv sysVar 	
	inert <- nlme::lme(sysVar ~ dist0:inert0 + dist0:inert1 + dist1:inert1 + dist1:inert0, random= ~ 1 | dyad, data= basedata, na.action=na.omit, method="ML", control=nlme::lmeControl(opt="optim"))
	names(inert$coefficients$fixed) <- c("intercept", aInert0name, pInert0name, aInert1name, pInert1name)
	
	# Coord indiv sysVar 
	coord <- nlme::lme(sysVar ~ dist0:coord0 + dist0:coord1 + dist1:coord1 + dist1:coord0, random= ~ 1 | dyad, data= basedata, na.action=na.omit, method="ML", control=nlme::lmeControl(opt="optim"))
	names(coord$coefficients$fixed) <- c("intercept", aCoord0name, pCoord0name, aCoord1name, pCoord1name)
	
	# InertCoord indiv sysVar 
	inertCoord <- nlme::lme(sysVar ~ dist0:inert0 + dist0:inert1 + dist1:inert1 + dist1:inert0 + dist0:coord0 + dist0:coord1 + dist1:coord1 + dist1:coord0, random= ~ 1 | dyad, data= basedata, na.action=na.omit, method="ML", control=nlme::lmeControl(opt="optim"))
	names(inertCoord$coefficients$fixed) <- c("intercept", aInert0name, pInert0name, aInert1name, pInert1name, aCoord0name, pCoord0name, aCoord1name, pCoord1name)

	basePred <- predict(base)
	inertPred <- predict(inert)
	coordPred <- predict(coord)
	inertCoordPred <- predict(inertCoord)
	
	obs <- basedata$sysVar
	baseR2 <- summary(lm(obs ~ basePred))$adj.r.squared
	inertR2 <- summary(lm(obs ~ inertPred))$adj.r.squared
	coordR2 <- summary(lm(obs ~ coordPred))$adj.r.squared
	inertCoordR2 <- summary(lm(obs ~ inertCoordPred))$adj.r.squared
	}

############### Plotting		

	min <- min(basedata$sysVar, na.rm=T)
	max <- max(basedata$sysVar, na.rm=T)
	
	ylabName <- paste(sysVarName, "observed", sep="_")
	xlabName <- paste(sysVarName, "predicted", sep="_")

	par(mfrow=c(2,1))
	hist(residuals(base))
	plot(basedata$sysVar ~ basePred, xlim=c(min, max), ylim=c(min, max), ylab=ylabName, xlab=xlabName, main="Baseline Model")
	
	par(mfrow=c(2,1))
	hist(residuals(inert))
	plot(basedata$sysVar ~ inertPred, xlim=c(min, max), ylim=c(min, max), ylab=ylabName, xlab=xlabName, main="Inertia Only Model")

	par(mfrow=c(2,1))
	hist(residuals(coord))
	plot(basedata$sysVar ~ coordPred, xlim=c(min, max), ylim=c(min, max), ylab=ylabName, xlab=xlabName, main="Coordination Only Model")
	
	par(mfrow=c(2,1))
	hist(residuals(inertCoord))
	plot(basedata$sysVar ~ inertCoordPred, xlim=c(min, max), ylim=c(min, max), ylab=ylabName, xlab=xlabName, main="Inertia-Coordination Model")
	
################ Collect results	
	models <- list(base=base, inert=inert, coord=coord, inertCoord=inertCoord)
	R2 <- list(baseR2=baseR2, inertR2=inertR2, coordR2=coordR2, inertCoordR2=inertCoordR2)
	
	output <- list(models=models, R2=R2)
}


#' Compares a baseline "intercepts only" model to one including the system variable for predicting the dynamic parameters of the inertia coordination model.
#' 
#' Multivariate correlated residuals dyadic models are used to predict the set of inertia-coordination parameters (inert0, inert1, coord0, coord1) from the system variable. The system variable can be either dyadic (sysVarType = "dyad"), where both partners have the same score (e.g., relationship length) or individual (sysVarType = "indiv"), where the partners can have different scores (e.g., age). If it is individual then both actor and partner effects of the system variable are included and are indicated by a_ and p_ preceding the variable names. 
#' 
#' @param basedata A dataframe containing the inertia-coordination parameter estimates produced by the "indivInertCoord" function.
#' @param sysVarType Whether the system variable is "dyad", which means both partners have the same socre, or "indiv" which means the partners can have different scores
#' @param dist0name A name for the level-0 of the distinguishing variable (e.g., "Women").
#' @param dist1name A name for the level-1 of the distinguishing variable (e.g., "Men").
#' 
#' @return The function returns a list including: 1) the gls objects containing the full results for each model (called "models"), and 2) adjusted R^2 information for each model (called "R2").  

#' @export
inertCoordSysVarIn <- function(basedata, sysVarType, dist0name, dist1name)
{
     basedata <- basedata[complete.cases(basedata), ] 
     
     # Names for intercepts
	inertInt0name <- paste("inert",dist0name, sep="_")
	inertInt1name <- paste("inert",dist1name, sep="_")
	coordInt0name <- paste("coord",dist0name, sep="_")
	coordInt1name <- paste("coord",dist1name, sep="_")
	
	# Names for individual model parameters
	aInert0name <- paste("a_Inert_sysVar",dist0name, sep="_")
	pInert0name <- paste("p_Inert_sysVar",dist0name, sep="_")
	aInert1name <- paste("a_Inert_sysVar",dist1name, sep="_")
	pInert1name <- paste("p_Inert_sysVar",dist1name, sep="_")
	aCoord0name <- paste("a_Coord_sysVar",dist0name, sep="_")
	aCoord1name <- paste("a_Coord_sysVar",dist1name, sep="_")
	pCoord0name <- paste("p_Coord_sysVar",dist0name, sep="_")
	pCoord1name <- paste("p_Coord_sysVar",dist1name, sep="_")

	# Names for dyad model parameters
	inertSysVar0name <- paste("inert_SysVar",dist0name, sep="_")
	inertSysVar1name <- paste("inert_SysVar",dist1name, sep="_")
	coordSysVar0name <- paste("coord_SysVar",dist0name, sep="_")
	coordSysVar1name <- paste("coord_SysVar",dist1name, sep="_")


      if(sysVarType != "indiv" & sysVarType != "dyad") 
   {
			stop("the sysVarType must be either indiv or dyad")
	}
 	else if(sysVarType=="indiv")
   {
	## Case with individual level sysVar

	data1 <- subset(basedata, select=c(dyad, sysVar, dist0))
	data2 <- stats::reshape(data1, idvar="dyad", timevar = "dist0", direction= "wide")
	data3 <- subset(basedata, select=c(dyad, inert0, inert1, coord0, coord1))
	data4 <- data3[!duplicated(data3$dyad), ]
	data5 <- plyr::join(data4, data2)
	colnames(data5) <- c("dyad", "paramEst1", "paramEst2", "paramEst3", "paramEst4", "sysVar1", "sysVar0")

	data6 <- stats::reshape(data5, varying=c("paramEst1", "paramEst2", "paramEst3", "paramEst4"), timevar="parameter", idvar="dyad", direction="long", sep="")
	data7 <- data6[complete.cases(data6), ] 

	data7$parameter <- factor(data7$parameter, levels = c(1,2,3,4), labels=c("inert0","inert1","coord0", "coord1"))
	data7$inert0 <- ifelse(data7$parameter == "inert0", 1, 0) # create indicator variables
	data7$inert1 <- ifelse(data7$parameter == "inert1", 1, 0) # create indicator variables
	data7$coord0 <- ifelse(data7$parameter == "coord0", 1, 0) # create indicator variables
	data7$coord1 <- ifelse(data7$parameter == "coord1", 1, 0) # create indicator variables

	base <- nlme::gls(paramEst ~ inert0 + inert1 + coord0 + coord1  -1,
	correlation=nlme::corSymm(form= ~ 1 |dyad), weights=nlme::varIdent(form = ~ 1|parameter), na.action=na.omit, method="ML", control=nlme::lmeControl(opt="optim"), data=data7)
	names(base$coefficients) <- c(inertInt0name, inertInt1name, coordInt0name, coordInt1name)
	

	sysVarIn <- nlme::gls(paramEst ~ inert0 + inert1 + coord0 + coord1 
				+ inert0:sysVar0 + inert0:sysVar1 + inert1:sysVar1 + inert1:sysVar0 + 
				coord0:sysVar0 + coord0:sysVar1 + coord1:sysVar1 + coord1:sysVar0 -1,
	correlation=nlme::corSymm(form= ~ 1 |dyad), weights=nlme::varIdent(form = ~ 1|parameter), na.action=na.omit, method="ML", control=nlme::lmeControl(opt="optim"), data=data7)
	names(sysVarIn$coefficients) <- c(inertInt0name, inertInt1name, coordInt0name, coordInt1name, aInert0name, pInert0name, aInert1name, pInert1name, aCoord0name, pCoord0name, aCoord1name, pCoord1name)

	
	obs <- data7$paramEst
			
	basePred <- predict(base)
	baseR2 <- summary(lm(obs ~ basePred))$adj.r.squared
	sysVarInPred <- predict(sysVarIn)
	sysVarInR2 <- summary(lm(obs ~ sysVarInPred))$adj.r.squared
	
	models <- list(base=base, sysVarIn=sysVarIn)
	R2 <- list(baseR2=baseR2, sysVarInR2=sysVarInR2)
				
	output <- list(models=models, R2=R2)
	}
	
	else if(sysVarType=="dyad")
	{
		## Case with dyadic level sysVar

	data1 <- basedata[!duplicated(basedata$dyad), ]
	data2 <- subset(data1, select=c(dyad, inert0, inert1, coord0, coord1, sysVar))
	colnames(data2) <- c("dyad","paramEst1", "paramEst2", "paramEst3", "paramEst4", "sysVar")
	data3 <- stats::reshape(data2, varying=c("paramEst1", "paramEst2", "paramEst3", "paramEst4"), timevar="parameter", idvar="dyad", direction="long", sep="")

	data3$parameter <- factor(data3$parameter, levels = c(1,2,3,4), labels=c("inert0","inert1","coord0", "coord1"))
	data3$inert0 <- ifelse(data3$parameter == "inert0", 1, 0) # create indicator variables
	data3$inert1 <- ifelse(data3$parameter == "inert1", 1, 0) # create indicator variables
	data3$coord0 <- ifelse(data3$parameter == "coord0", 1, 0) # create indicator variables
	data3$coord1 <- ifelse(data3$parameter == "coord1", 1, 0) # create indicator variables

	base <- nlme::gls(paramEst ~ inert0 + inert1 + coord0 + coord1  -1, correlation=nlme::corSymm(form= ~ 1 |dyad), weights=nlme::varIdent(form = ~ 1|parameter), na.action=na.omit, method="ML", control=nlme::lmeControl(opt="optim"), data=data3)
	names(base$coefficients) <- c(inertInt0name, inertInt1name, coordInt0name, coordInt1name)


	sysVarIn <- nlme::gls(paramEst ~ inert0 + inert1 + coord0 + coord1 
				+ sysVar:inert0 + sysVar:inert1 + sysVar:coord0 + sysVar:coord1 -1, correlation=nlme::corSymm(form= ~ 1 |dyad), weights=nlme::varIdent(form = ~ 1|parameter), na.action=na.omit, method="ML", control=nlme::lmeControl(opt="optim"), data=data3)
	names(sysVarIn$coefficients) <- c(inertInt0name, inertInt1name, coordInt0name, coordInt1name, inertSysVar0name, inertSysVar1name, coordSysVar0name, coordSysVar1name)
	
	obs <- data3$paramEst
			
	basePred <- predict(base)
	baseR2 <- summary(lm(obs ~ basePred))$adj.r.squared
	sysVarInPred <- predict(sysVarIn)
	sysVarInR2 <- summary(lm(obs ~ sysVarInPred))$adj.r.squared
	
	models <- list(base=base, sysVarIn=sysVarIn)
	R2 <- list(baseR2=baseR2, sysVarInR2=sysVarInR2)
				
	output <- list(models=models, R2=R2)
	} 	
}
















