

######## This file includes all the functions needed for an inertia-coordination analysis

#' Estimates an "inertia only" model for each dyad.
#' 
#' The observed state variables (with linear trends removed) are predicted from each person's intercept and each person's own observed state variable lagged at the amount specified during the dataPrep step (again with linear trends removed)
#'
#' @param basedata A dataframe that was produced with the "dataPrep" function.
#' @param dist0name A name for the 0-level of the distinguishing variable (e.g., "Women").
#' @param dist1name A name for the 1-level of the distinguishing variable (e.g., "Men").
#' @param obsName A name for the observed state variables being plotted (e.g., "Emotional Experience").
#' 
#' @return The function returns a list including: 1) the adjusted R^2 for the model for each dyad (called "r2"), 2) the parameter estimates for the model for each dyad (called "paramData", for use in either predicting, or being predicted by, the system variable), and 3) plots of the predicted values against the observed values for each dyad (called "plots"). The plots are also written to the working directory as a pdf file called "inertPlots.pdf"

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

	results <- list(r2=r2, paramData=paramData, plots=plots)
}


#' Estimates a "coordination only" model for each dyad.
#' 
#' The observed state variables (with linear trends removed) are predicted from each person's intercept and each person's partner's concurrent state variable (again with linear trends removed).
#' 
#' @param basedata A dataframe that was produced with the "dataPrep" function.
#' @param dist0name A name for the 0-level of the distinguishing variable (e.g., "Women").
#' @param dist1name A name for the 1-level of the distinguishing variable (e.g., "Men").
#' @param obsName A name for the observed state variables being plotted (e.g., "Emotional Experience").
#' 
#' @return The function returns a list including: 1) the adjusted R^2 for the model for each dyad (called "r2", 2) the parameter estimates for the model for each dyad (called "paramData", for use in either predicting, or being predicted by, the system variable), and 3) plots of the predicted values against the observed values for each dyad (called "plots"). The plots are also written to the working directory as a pdf file called "coordPlots.pdf"

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

	results <- list(r2=r2, paramData=paramData, plots=plots)
}


#' Estimates the full "inertia-coordination" model for each dyad.
#' 
#' The observed state variables (with linear trends removed) are predicted from each person's intercept, each person's own observed state variable lagged at the amount specified during the dataPrep step (again with linear trends removed), and each person's partner's concurrent state variable (again with linear trends removed).
#' 
#' @param basedata A dataframe that was produced with the "dataPrep" function.
#' @param dist0name A name for the 0-level of the distinguishing variable (e.g., "Women").
#' @param dist1name A name for the 1-level of the distinguishing variable (e.g., "Men").
#' @param obsName A name for the observed state variables being plotted (e.g., "Emotional Experience").
#' 
#' @return The function returns a list including: 1) the adjusted R^2 for the model for each dyad (called "r2", 2) the parameter estimates for the model for each dyad (called "paramData", for use in either predicting, or being predicted by, the system variable), and 3) plots of the predicted values against the observed values for each dyad (called "plots"). The plots are also written to the working directory as a pdf file called "inertCoordPlots.pdf"

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

	results <- list(r2=r2, paramData=paramData, plots=plots)
}


#' Compares inertia-coordination models for predicting the system variable from the dynamic parameters.
#' 
#' The dynamic parameters used in these models come from a set including both people's inertia (inert0 and inert1) and coordination (coord0 and coord1) estimates. The 3 models compared are the inertia-only (inert0 + inert1), coordination-only (coord0 + coord1) and the full inertia-coordination (inert0 + inert1 + coord0 + coord1) models. These are estimated separately for each level of the distinguishing variable. 
#' 
#' @param basedata A dataframe containing the inertia-coordination parameter estimates produced by the "indivInertCoord" function.
#' @param dist0name A name for the 0-level of the distinguishing variable (e.g., "Women").
#' @param dist1name A name for the 1-level of the distinguishing variable (e.g., "Men").
#' @param sysVarName A name for the system variable being predicted (e.g., "Satisfaction").
#' 
#' @return The function returns a list including: 1) the lm objects containing the full results for each model for each level of the distinguishing variable (called "models0" and "models1"), 2) anova output for each model for each level of the distinguishing variable (called "anovas0" and "anovas1"), 3) summary output for each model for each level of the distinguishing variable (called "summaries0" and "summaries1") and 4) adjusted R^2 information for each model for each level of the distinguishing variable (called "adjustR20" and "adjustR21"). The function also displays histograms of the residuals and plots of the predicted values against observed values for each model. 

#' @export
inertCoordSysVarOutCompare <- function(basedata, dist0name, dist1name, sysVarName)
{
	basedata <- basedata[complete.cases(basedata), ] # remove NAs in order to compare nested models	
	basedata0 <- subset(basedata, dist0==1)
	basedata1 <- subset(basedata, dist0==0)

	inert0name <- paste("inert",dist0name, sep="_")
	inert1name <- paste("inert",dist1name, sep="_")
	coord0name <- paste("coord",dist0name, sep="_")
	coord1name <- paste("coord",dist1name, sep="_")

	inert0 <- lm(sysVar ~ inert0 + inert1 , data= basedata0)
	names(inert0$coefficients) <- c("intercept", inert0name, inert1name)
	inert1 <- lm(sysVar ~ inert0 + inert1 , data= basedata1)
	names(inert1$coefficients) <- c("intercept", inert0name, inert1name)

	coord0 <- lm(sysVar ~ coord0 + coord1 , data= basedata0)
	names(coord0$coefficients) <- c("intercept", coord0name, coord1name)
	coord1 <- lm(sysVar ~ coord0 + coord1 , data= basedata1)
	names(coord1$coefficients) <- c("intercept", coord0name, coord1name)
	
	inertCoord0 <- lm(sysVar ~ inert0 + coord0 + inert1 + coord1 , data= basedata0)
	names(inertCoord0$coefficients) <- c("intercept", inert0name, coord0name, inert1name, coord1name)
	inertCoord1 <- lm(sysVar ~ inert0 + coord0 + inert1 + coord1 , data= basedata1)
	names(inertCoord1$coefficients) <- c("intercept", inert0name, coord0name, inert1name, coord1name)

	inert0Anova <- car::Anova(inert0, type="III")
	inert0Summary <- summary(inert0)
	coord0Anova <- car::Anova(coord0, type="III")
	coord0Summary <- summary(coord0)
	inertCoord0Anova <- car::Anova(inertCoord0, type="III")
	inertCoord0Summary <- summary(inertCoord0)
	
	inert1Anova <- car::Anova(inert1, type="III")
	inert1Summary <- summary(inert1)
	coord1Anova <- car::Anova(coord1, type="III")
	coord1Summary <- summary(coord1)
	inertCoord1Anova <- car::Anova(inertCoord1, type="III")
	inertCoord1Summary <- summary(inertCoord1)	
	
	inert0Pred <- predict(inert0)
	inert0R2 <- summary(inert0)$adj.r.squared
	coord0Pred <- predict(coord0)
	coord0R2 <- summary(coord0)$adj.r.squared
	inertCoord0Pred <- predict(inertCoord0)
	inertCoord0R2 <- summary(inertCoord0)$adj.r.squared
	
	inert1Pred <- predict(inert1)
	inert1R2 <- summary(inert1)$adj.r.squared
	coord1Pred <- predict(coord1)
	coord1R2 <- summary(coord1)$adj.r.squared
	inertCoord1Pred <- predict(inertCoord1)
	inertCoord1R2 <- summary(inertCoord1)$adj.r.squared	

	min0 <- min(basedata0$sysVar, na.rm=T)
	max0 <- max(basedata0$sysVar, na.rm=T)
	min1 <- min(basedata1$sysVar, na.rm=T)
	max1 <- max(basedata1$sysVar, na.rm=T)
	
	ylabName <- paste(sysVarName, "obs", sep="_")
	xlabName <- paste(sysVarName, "pred", sep="_")
	inert0name <- paste("Inert", dist0name, sep="_")
	inert1name <- paste("Inert", dist1name, sep="_")
	coord0name <- paste("Coord", dist0name, sep="_")
	coord1name <- paste("Coord", dist1name, sep="_")
	inertCoord0name <- paste("InertCoord", dist0name, sep="_")
	inertCoord1name <- paste("InertCoord", dist1name, sep="_")

	par(mfrow=c(2,1))
	hist(residuals(inert0))
	plot(basedata0$sysVar ~ inert0Pred, xlim=c(min0, max0), ylim=c(min0, max0), ylab=ylabName, xlab=xlabName, main=inert0name)
	par(mfrow=c(2,1))
	hist(residuals(coord0))
	plot(basedata0$sysVar ~ coord0Pred, xlim=c(min0, max0), ylim=c(min0, max0), ylab=ylabName, xlab=xlabName, main=coord0name)
	par(mfrow=c(2,1))
	hist(residuals(inertCoord0))
	plot(basedata0$sysVar ~ inertCoord0Pred, xlim=c(min0, max0), ylim=c(min0, max0), ylab=ylabName, xlab=xlabName, main=inertCoord0name)
	
	par(mfrow=c(2,1))
	hist(residuals(inert1))
	plot(basedata1$sysVar ~ inert1Pred, xlim=c(min1, max1), ylim=c(min1, max1), ylab=ylabName, xlab=xlabName, main=inert1name)
	par(mfrow=c(2,1))
	hist(residuals(coord1))
	plot(basedata1$sysVar ~ coord1Pred, xlim=c(min1, max1), ylim=c(min1, max1), ylab=ylabName, xlab=xlabName, main=coord1name)
	par(mfrow=c(2,1))
	hist(residuals(inertCoord1))
	plot(basedata1$sysVar ~ inertCoord1Pred, xlim=c(min1, max1), ylim=c(min1, max1), ylab=ylabName, xlab=xlabName, main=inertCoord1name)
	
	
	models0 <- list(inert0=inert0, coord0=coord0, inertCoord0=inertCoord0)
	models1 <- list(inert1=inert1, coord1=coord1, inertCoord1=inertCoord1)
	anovas0 <- list(inert0Anova=inert0Anova, coord0Anova=coord0Anova, inertCoord0Anova=inertCoord0Anova)
	anovas1 <- list(inert1Anova=inert1Anova, coord1Anova=coord1Anova, inertCoord1Anova=inertCoord1Anova)
	summaries0 <- list(inert0Summary=inert0Summary, coord0Summary=coord0Summary, inertCoord0Summary=inertCoord0Summary)
	summaries1 <- list(inert1Summary=inert1Summary, coord1Summary=coord1Summary, inertCoord1Summary=inertCoord1Summary)
	r20 <- list(inert0R2=inert0R2, coord0R2=coord0R2, inertCoord0R2=inertCoord0R2)
	r21 <- list(inert1R2=inert1R2, coord1R2=coord1R2, inertCoord1R2=inertCoord1R2)
	
	output <- list(models0=models0, models1=models1, anovas0=anovas0, summaries0=summaries0, adjustR20=r20, anovas1=anovas1, summaries1=summaries1, adjustR21=r21)
}


#' Plots of the system variable predicted from the inertia parameter estimates.
#'
#' Displays bar plots for one level of the distinguishing variable, showing model predicted means and standard errors of the system variable at low and high levels of a person's own (actor) and partner's (partner) inertia parameter estimates.
#' 
#' @param basedata A dataframe containing the inertia-coordination parameter estimates produced by the "indivInertCoord" function.
#' @param centInert0 A vector of low, medium and high centering values for the inertia parameter estimates for the 0-level of the distinguishing variable (e.g., values that indicate typical low, medium and high values of inertia for the 0-level partner).
#' @param centInert1 A vector of low, medium and high centering values for the inertia parameter estimates for the 1-level of the distinguishing variable (e.g., values that indicate typical low, medium and high values of inertia for the 1-level partner).
#' @param sysVarName A name for the system variable being predicted (e.g., "Satisfaction").
#' @param dist A number indicating which level of the distinguishing variable to produce plots for. Must be either 0 or 1.
#' @param dist0name A name for the 0-level of the distinguishing variable (e.g., "Women").
#' @param dist1name A name for the 1-level of the distinguishing variable (e.g., "Men").

#' @import ggplot2
#' @export
inertSysVarOutPlots <- function(basedata, centInert0, centInert1, sysVarName, dist, dist0name, dist1name)
{
	if(dist==0)
	{
		basedata <- basedata[basedata$dist0==1, ]
		distName <- dist0name
	} else if (dist==1)
		{
			basedata <- basedata[basedata$dist0 == 0, ]
			distName <- dist1name
		}else {cat("\n error: dist must be 0 or 1\n")
			 	stop(call.=F)}
	
	if(dist==0)
	{
	inertAData <- basedata
	inertAData$inert0L <- basedata$inert0 - centInert0[1]
	inertAData$inert1M <- basedata$inert1 - centInert1[2]
	inertAData$inert0H <- basedata$inert0 - centInert0[3]
	inertAData$inert1M <- basedata$inert1 - centInert1[2]
	
	inertPData <- basedata
	inertPData$inert0M <- basedata$inert0 - centInert0[2]
	inertPData$inert1L <- basedata$inert1 - centInert1[1]
	inertPData$inert0M <- basedata$inert0 - centInert0[2]
	inertPData$inert1H <- basedata$inert1 - centInert1[3]
	
	inertAL <- lm(sysVar ~ inert0L + inert1M, data=inertAData)
	inertAH <- lm(sysVar ~ inert0H + inert1M, data=inertAData)
	inertPL <- lm(sysVar ~ inert0M + inert1L, data=inertPData)
	inertPH <- lm(sysVar ~ inert0M + inert1H, data=inertPData)
	
	estInertAL <- summary(inertAL)$coefficients[1,1]
	errInertAL <- summary(inertAL)$coefficients[1,2]
	estInertAH <- summary(inertAH)$coefficients[1,1]
	errInertAH <- summary(inertAH)$coefficients[1,2]
	
	estInertPL <- summary(inertPL)$coefficients[1,1]
	errInertPL <- summary(inertPL)$coefficients[1,2]
	estInertPH <- summary(inertPH)$coefficients[1,1]
	errInertPH <- summary(inertPH)$coefficients[1,2]

	inertALdata <- c(1,1,estInertAL,errInertAL)
	inertAHdata <- c(1,2,estInertAH,errInertAH)
	inertPLdata <- c(2,1,estInertPL,errInertPL)
	inertPHdata <- c(2,2,estInertPH,errInertPH)
		
	plotData <- as.data.frame(rbind(inertALdata, inertAHdata, inertPLdata, inertPHdata))
	names(plotData) <- c("Parameter","Level","sysVar","SE")
	plotData$Parameter <- factor(plotData$Parameter, labels=c("Actor", "Partner"))
	plotData$Level <- factor(plotData$Level, labels=c("Low","High"))
	plotData$errMin <- plotData[ ,3] - plotData[ ,4]
	plotData$errMax <- plotData[ ,3] + plotData[ ,4]
		
	Ymin <- min(basedata$sysVar, na.rm=T)
	Ymax <- max(basedata$sysVar, na.rm=T)
			
	inert <- ggplot(plotData, aes(x=Parameter, y=sysVar, fill=Level)) +
		geom_bar(position=position_dodge(), stat="identity") +
		geom_errorbar(aes(ymin=errMin, ymax=errMax), width=.2, position=position_dodge(.9)) +
		coord_cartesian(ylim=c(Ymin, Ymax)) +
		labs(title="Inertia", subtitle= distName) +
		ylab(sysVarName) +
		theme(plot.title=element_text(size=12))
	
	} else if (dist==1)
	
	{
	inertAData <- basedata
	inertAData$inert0M <- basedata$inert0 - centInert0[2]
	inertAData$inert1L <- basedata$inert1 - centInert1[1]
	inertAData$inert0M <- basedata$inert0 - centInert0[2]
	inertAData$inert1H <- basedata$inert1 - centInert1[3]
	
	inertPData <- basedata
	inertPData$inert0L <- basedata$inert0 - centInert0[1]
	inertPData$inert1M <- basedata$inert1 - centInert1[2]
	inertPData$inert0H <- basedata$inert0 - centInert0[3]
	inertPData$inert1M <- basedata$inert1 - centInert1[2]
	
	inertAL <- lm(sysVar ~ inert0M + inert1L, data=inertAData)
	inertAH <- lm(sysVar ~ inert0M + inert1H, data=inertAData)
	inertPL <- lm(sysVar ~ inert0L + inert1M, data=inertPData)
	inertPH <- lm(sysVar ~ inert0H + inert1M, data=inertPData)
	
	estInertAL <- summary(inertAL)$coefficients[1,1]
	errInertAL <- summary(inertAL)$coefficients[1,2]
	estInertAH <- summary(inertAH)$coefficients[1,1]
	errInertAH <- summary(inertAH)$coefficients[1,2]
	
	estInertPL <- summary(inertPL)$coefficients[1,1]
	errInertPL <- summary(inertPL)$coefficients[1,2]
	estInertPH <- summary(inertPH)$coefficients[1,1]
	errInertPH <- summary(inertPH)$coefficients[1,2]

	inertALdata <- c(1,1,estInertAL,errInertAL)
	inertAHdata <- c(1,2,estInertAH,errInertAH)
	inertPLdata <- c(2,1,estInertPL,errInertPL)
	inertPHdata <- c(2,2,estInertPH,errInertPH)
		
	plotData <- as.data.frame(rbind(inertALdata, inertAHdata, inertPLdata, inertPHdata))
	names(plotData) <- c("Parameter","Level","sysVar","SE")
	plotData$Parameter <- factor(plotData$Parameter, labels=c("Actor", "Partner"))
	plotData$Level <- factor(plotData$Level, labels=c("Low","High"))
	plotData$errMin <- plotData[ ,3] - plotData[ ,4]
	plotData$errMax <- plotData[ ,3] + plotData[ ,4]
		
	Ymin <- min(basedata$sysVar, na.rm=T)
	Ymax <- max(basedata$sysVar, na.rm=T)
			
	inert <- ggplot(plotData, aes(x=Parameter, y=sysVar, fill=Level)) +
		geom_bar(position=position_dodge(), stat="identity") +
		geom_errorbar(aes(ymin=errMin, ymax=errMax), width=.2, position=position_dodge(.9)) +
		coord_cartesian(ylim=c(Ymin, Ymax)) +
		labs(title="Inertia", subtitle= distName) +
		ylab(sysVarName) +
		theme(plot.title=element_text(size=12))	
	}
	plots <- inert
}

#' Plots of the system variable predicted from the coordination parameter estimates.
#'
#' Displays bar plots for one level of the distinguishing variable, showing model predicted means and standard errors of the system variable at low and high levels of a person's own (actor) and partner's (partner) coordination parameter estimates.
#' 
#' @param basedata A dataframe containing the inertia-coordination parameter estimates produced by the "indivInertCoord" function.
#' @param centCoord0 A vector of low, medium and high centering values for the coordination parameter estimates for the 0-level of the distinguishing variable (e.g., values that indicate typical low, medium and high values of coordination for the 0-level partner).
#' @param centCoord1 A vector of low, medium and high centering values for the coordination parameter estimates for the 1-level of the distinguishing variable (e.g., values that indicate typical low, medium and high values of coordination for the 1-level partner).
#' @param sysVarName A name for the system variable being predicted (e.g., "Satisfaction").
#' @param dist A number indicating which level of the distinguishing variable to produce plots for. Must be either 0 or 1.
#' @param dist0name A name for the 0-level of the distinguishing variable (e.g., "Women").
#' @param dist1name A name for the 1-level of the distinguishing variable (e.g., "Men").

#' @import ggplot2
#' @export
coordSysVarOutPlots <- function(basedata, centCoord0, centCoord1, sysVarName, dist, dist0name, dist1name)
{
	if(dist==0)
	{
		basedata <- basedata[basedata$dist0==1, ]
		distName <- dist0name
	} else if (dist==1)
		{
			basedata <- basedata[basedata$dist0 == 0, ]
			distName <- dist1name
		}else {cat("\n error: dist must be 0 or 1\n")
			 	stop(call.=F)}
	
	if(dist==0)
	{
	coordAData <- basedata
	coordAData$coord0L <- basedata$coord0 - centCoord0[1]
	coordAData$coord1M <- basedata$coord1 - centCoord1[2]
	coordAData$coord0H <- basedata$coord0 - centCoord0[3]
	coordAData$coord1M <- basedata$coord1 - centCoord1[2]
	
	coordPData <- basedata
	coordPData$coord0M <- basedata$coord0 - centCoord0[2]
	coordPData$coord1L <- basedata$coord1 - centCoord1[1]
	coordPData$coord0M <- basedata$coord0 - centCoord0[2]
	coordPData$coord1H <- basedata$coord1 - centCoord1[3]
	
	coordAL <- lm(sysVar ~ coord0L + coord1M, data=coordAData)
	coordAH <- lm(sysVar ~ coord0H + coord1M, data=coordAData)
	coordPL <- lm(sysVar ~ coord0M + coord1L, data=coordPData)
	coordPH <- lm(sysVar ~ coord0M + coord1H, data=coordPData)
	
	estcoordAL <- summary(coordAL)$coefficients[1,1]
	errcoordAL <- summary(coordAL)$coefficients[1,2]
	estcoordAH <- summary(coordAH)$coefficients[1,1]
	errcoordAH <- summary(coordAH)$coefficients[1,2]
	
	estcoordPL <- summary(coordPL)$coefficients[1,1]
	errcoordPL <- summary(coordPL)$coefficients[1,2]
	estcoordPH <- summary(coordPH)$coefficients[1,1]
	errcoordPH <- summary(coordPH)$coefficients[1,2]

	coordALdata <- c(1,1,estcoordAL,errcoordAL)
	coordAHdata <- c(1,2,estcoordAH,errcoordAH)
	coordPLdata <- c(2,1,estcoordPL,errcoordPL)
	coordPHdata <- c(2,2,estcoordPH,errcoordPH)
		
	plotData <- as.data.frame(rbind(coordALdata, coordAHdata, coordPLdata, coordPHdata))
	names(plotData) <- c("Parameter","Level","sysVar","SE")
	plotData$Parameter <- factor(plotData$Parameter, labels=c("Actor", "Partner"))
	plotData$Level <- factor(plotData$Level, labels=c("Low","High"))
	plotData$errMin <- plotData[ ,3] - plotData[ ,4]
	plotData$errMax <- plotData[ ,3] + plotData[ ,4]
		
	Ymin <- min(basedata$sysVar, na.rm=T)
	Ymax <- max(basedata$sysVar, na.rm=T)
			
	coord <- ggplot(plotData, aes(x=Parameter, y=sysVar, fill=Level)) +
		geom_bar(position=position_dodge(), stat="identity") +
		geom_errorbar(aes(ymin=errMin, ymax=errMax), width=.2, position=position_dodge(.9)) +
		coord_cartesian(ylim=c(Ymin, Ymax)) +
		labs(title="Coordination", subtitle= distName) +
		ylab(sysVarName) +
		theme(plot.title=element_text(size=12))
	
	} else if (dist==1)
	
	{
	coordAData <- basedata
	coordAData$coord0M <- basedata$coord0 - centCoord0[2]
	coordAData$coord1L <- basedata$coord1 - centCoord1[1]
	coordAData$coord0M <- basedata$coord0 - centCoord0[2]
	coordAData$coord1H <- basedata$coord1 - centCoord1[3]
	
	coordPData <- basedata
	coordPData$coord0L <- basedata$coord0 - centCoord0[1]
	coordPData$coord1M <- basedata$coord1 - centCoord1[2]
	coordPData$coord0H <- basedata$coord0 - centCoord0[3]
	coordPData$coord1M <- basedata$coord1 - centCoord1[2]
	
	coordAL <- lm(sysVar ~ coord0M + coord1L, data=coordAData)
	coordAH <- lm(sysVar ~ coord0M + coord1H, data=coordAData)
	coordPL <- lm(sysVar ~ coord0L + coord1M, data=coordPData)
	coordPH <- lm(sysVar ~ coord0H + coord1M, data=coordPData)
	
	estcoordAL <- summary(coordAL)$coefficients[1,1]
	errcoordAL <- summary(coordAL)$coefficients[1,2]
	estcoordAH <- summary(coordAH)$coefficients[1,1]
	errcoordAH <- summary(coordAH)$coefficients[1,2]
	
	estcoordPL <- summary(coordPL)$coefficients[1,1]
	errcoordPL <- summary(coordPL)$coefficients[1,2]
	estcoordPH <- summary(coordPH)$coefficients[1,1]
	errcoordPH <- summary(coordPH)$coefficients[1,2]

	coordALdata <- c(1,1,estcoordAL,errcoordAL)
	coordAHdata <- c(1,2,estcoordAH,errcoordAH)
	coordPLdata <- c(2,1,estcoordPL,errcoordPL)
	coordPHdata <- c(2,2,estcoordPH,errcoordPH)
		
	plotData <- as.data.frame(rbind(coordALdata, coordAHdata, coordPLdata, coordPHdata))
	names(plotData) <- c("Parameter","Level","sysVar","SE")
	plotData$Parameter <- factor(plotData$Parameter, labels=c("Actor", "Partner"))
	plotData$Level <- factor(plotData$Level, labels=c("Low","High"))
	plotData$errMin <- plotData[ ,3] - plotData[ ,4]
	plotData$errMax <- plotData[ ,3] + plotData[ ,4]
		
	Ymin <- min(basedata$sysVar, na.rm=T)
	Ymax <- max(basedata$sysVar, na.rm=T)
			
	coord <- ggplot(plotData, aes(x=Parameter, y=sysVar, fill=Level)) +
		geom_bar(position=position_dodge(), stat="identity") +
		geom_errorbar(aes(ymin=errMin, ymax=errMax), width=.2, position=position_dodge(.9)) +
		coord_cartesian(ylim=c(Ymin, Ymax)) +
		labs(title="Coordination", subtitle= distName) +
		ylab(sysVarName) +
		theme(plot.title=element_text(size=12))
	
	}
	plots <- coord
}

#' Plots of the system variable predicted from the inertia and coordination parameter estimates.
#'
#' Displays bar plots for one level of the distinguishing variable, showing model predicted means and standard errors of the system variable at low and high levels of a person's own (actor) and partner's (partner) inertia and coordination parameter estimates.
#' 
#' @param basedata A dataframe containing the inertia-coordination parameter estimates produced by the "indivInertCoord" function.
#' @param centInert0 A vector of low, medium and high centering values for the inertia parameter estimates for the 0-level of the distinguishing variable (e.g., values that indicate typical low, medium and high values of inertia for the 0-level partner).
#' @param centInert1 A vector of low, medium and high centering values for the inertia parameter estimates for the 1-level of the distinguishing variable (e.g., values that indicate typical low, medium and high values of inertia for the 1-level partner).
#' @param centCoord0 A vector of low, medium and high centering values for the coordination parameter estimates for the 0-level of the distinguishing variable (e.g., values that indicate typical low, medium and high values of coordination for the 0-level partner).
#' @param centCoord1 A vector of low, medium and high centering values for the coordination parameter estimates for the 1-level of the distinguishing variable (e.g., values that indicate typical low, medium and high values of coordination for the 1-level partner).
#' @param sysVarName A name for the system variable being predicted (e.g., "Satisfaction").
#' @param dist A number indicating which level of the distinguishing variable to produce plots for. Must be either 0 or 1.
#' @param dist0name A name for the 0-level of the distinguishing variable (e.g., "Women").
#' @param dist1name A name for the 1-level of the distinguishing variable (e.g., "Men").

#' @import ggplot2
#' @export

inertCoordSysVarOutPlots <- function(basedata, centInert0, centInert1, centCoord0, centCoord1, sysVarName, dist, dist0name, dist1name)
{
	if(dist==0)
	{
		basedata <- basedata[basedata$dist0==1, ]
		distName <- dist0name
	} else if (dist==1)
		{
			basedata <- basedata[basedata$dist0 == 0, ]
			distName <- dist1name
		}else {cat("\n error: dist must be 0 or 1\n")
			 	stop(call.=F)}
	
	if(dist==0)
	{
	
	inertAData <- basedata
	inertAData$inert0L <- basedata$inert0 - centInert0[1]
	inertAData$inert1M <- basedata$inert1 - centInert1[2]
	inertAData$inert0H <- basedata$inert0 - centInert0[3]
	inertAData$inert1M <- basedata$inert1 - centInert1[2]
	
	inertPData <- basedata
	inertPData$inert0M <- basedata$inert0 - centInert0[2]
	inertPData$inert1L <- basedata$inert1 - centInert1[1]
	inertPData$inert0M <- basedata$inert0 - centInert0[2]
	inertPData$inert1H <- basedata$inert1 - centInert1[3]
	
	inertAL <- lm(sysVar ~ inert0L + inert1M, data=inertAData)
	inertAH <- lm(sysVar ~ inert0H + inert1M, data=inertAData)
	inertPL <- lm(sysVar ~ inert0M + inert1L, data=inertPData)
	inertPH <- lm(sysVar ~ inert0M + inert1H, data=inertPData)
	
	estInertAL <- summary(inertAL)$coefficients[1,1]
	errInertAL <- summary(inertAL)$coefficients[1,2]
	estInertAH <- summary(inertAH)$coefficients[1,1]
	errInertAH <- summary(inertAH)$coefficients[1,2]
	
	estInertPL <- summary(inertPL)$coefficients[1,1]
	errInertPL <- summary(inertPL)$coefficients[1,2]
	estInertPH <- summary(inertPH)$coefficients[1,1]
	errInertPH <- summary(inertPH)$coefficients[1,2]

	inertALdata <- c(1,1,estInertAL,errInertAL)
	inertAHdata <- c(1,2,estInertAH,errInertAH)
	inertPLdata <- c(2,1,estInertPL,errInertPL)
	inertPHdata <- c(2,2,estInertPH,errInertPH)
		
	plotDataInert <- as.data.frame(rbind(inertALdata, inertAHdata, inertPLdata, inertPHdata))
	names(plotDataInert) <- c("Parameter","Level","sysVar","SE")
	plotDataInert$Parameter <- factor(plotDataInert$Parameter, labels=c("Actor", "Partner"))
	plotDataInert$Level <- factor(plotDataInert$Level, labels=c("Low","High"))
	plotDataInert$errMin <- plotDataInert[ ,3] - plotDataInert[ ,4]
	plotDataInert$errMax <- plotDataInert[ ,3] + plotDataInert[ ,4]
		
	Ymin <- min(basedata$sysVar, na.rm=T)
	Ymax <- max(basedata$sysVar, na.rm=T)
			
	inert <- ggplot(plotDataInert, aes(x=Parameter, y=sysVar, fill=Level)) +
		geom_bar(position=position_dodge(), stat="identity") +
		geom_errorbar(aes(ymin=errMin, ymax=errMax), width=.2, position=position_dodge(.9)) +
		coord_cartesian(ylim=c(Ymin, Ymax)) +
		labs(title="Inertia", subtitle= distName) +
		ylab(sysVarName) +
		theme(plot.title=element_text(size=12))
	

	coordAData <- basedata
	coordAData$coord0L <- basedata$coord0 - centCoord0[1]
	coordAData$coord1M <- basedata$coord1 - centCoord1[2]
	coordAData$coord0H <- basedata$coord0 - centCoord0[3]
	coordAData$coord1M <- basedata$coord1 - centCoord1[2]
	
	coordPData <- basedata
	coordPData$coord0M <- basedata$coord0 - centCoord0[2]
	coordPData$coord1L <- basedata$coord1 - centCoord1[1]
	coordPData$coord0M <- basedata$coord0 - centCoord0[2]
	coordPData$coord1H <- basedata$coord1 - centCoord1[3]
	
	coordAL <- lm(sysVar ~ coord0L + coord1M, data=coordAData)
	coordAH <- lm(sysVar ~ coord0H + coord1M, data=coordAData)
	coordPL <- lm(sysVar ~ coord0M + coord1L, data=coordPData)
	coordPH <- lm(sysVar ~ coord0M + coord1H, data=coordPData)
	
	estcoordAL <- summary(coordAL)$coefficients[1,1]
	errcoordAL <- summary(coordAL)$coefficients[1,2]
	estcoordAH <- summary(coordAH)$coefficients[1,1]
	errcoordAH <- summary(coordAH)$coefficients[1,2]
	
	estcoordPL <- summary(coordPL)$coefficients[1,1]
	errcoordPL <- summary(coordPL)$coefficients[1,2]
	estcoordPH <- summary(coordPH)$coefficients[1,1]
	errcoordPH <- summary(coordPH)$coefficients[1,2]

	coordALdata <- c(1,1,estcoordAL,errcoordAL)
	coordAHdata <- c(1,2,estcoordAH,errcoordAH)
	coordPLdata <- c(2,1,estcoordPL,errcoordPL)
	coordPHdata <- c(2,2,estcoordPH,errcoordPH)
		
	plotDataCoord <- as.data.frame(rbind(coordALdata, coordAHdata, coordPLdata, coordPHdata))
	names(plotDataCoord) <- c("Parameter","Level","sysVar","SE")
	plotDataCoord$Parameter <- factor(plotDataCoord$Parameter, labels=c("Actor", "Partner"))
	plotDataCoord$Level <- factor(plotDataCoord$Level, labels=c("Low","High"))
	plotDataCoord$errMin <- plotDataCoord[ ,3] - plotDataCoord[ ,4]
	plotDataCoord$errMax <- plotDataCoord[ ,3] + plotDataCoord[ ,4]
		
	Ymin <- min(basedata$sysVar, na.rm=T)
	Ymax <- max(basedata$sysVar, na.rm=T)
			
	coord <- ggplot(plotDataCoord, aes(x=Parameter, y=sysVar, fill=Level)) +
		geom_bar(position=position_dodge(), stat="identity") +
		geom_errorbar(aes(ymin=errMin, ymax=errMax), width=.2, position=position_dodge(.9)) +
		coord_cartesian(ylim=c(Ymin, Ymax)) +
		labs(title="Coordination", subtitle= distName) +
		ylab(sysVarName) +
		theme(plot.title=element_text(size=12))
	
	} else if (dist==1)
	
	{
	
	inertAData <- basedata
	inertAData$inert0M <- basedata$inert0 - centInert0[2]
	inertAData$inert1L <- basedata$inert1 - centInert1[1]
	inertAData$inert0M <- basedata$inert0 - centInert0[2]
	inertAData$inert1H <- basedata$inert1 - centInert1[3]
	
	inertPData <- basedata
	inertPData$inert0L <- basedata$inert0 - centInert0[1]
	inertPData$inert1M <- basedata$inert1 - centInert1[2]
	inertPData$inert0H <- basedata$inert0 - centInert0[3]
	inertPData$inert1M <- basedata$inert1 - centInert1[2]
	
	inertAL <- lm(sysVar ~ inert0M + inert1L, data=inertAData)
	inertAH <- lm(sysVar ~ inert0M + inert1H, data=inertAData)
	inertPL <- lm(sysVar ~ inert0L + inert1M, data=inertPData)
	inertPH <- lm(sysVar ~ inert0H + inert1M, data=inertPData)
	
	estInertAL <- summary(inertAL)$coefficients[1,1]
	errInertAL <- summary(inertAL)$coefficients[1,2]
	estInertAH <- summary(inertAH)$coefficients[1,1]
	errInertAH <- summary(inertAH)$coefficients[1,2]
	
	estInertPL <- summary(inertPL)$coefficients[1,1]
	errInertPL <- summary(inertPL)$coefficients[1,2]
	estInertPH <- summary(inertPH)$coefficients[1,1]
	errInertPH <- summary(inertPH)$coefficients[1,2]

	inertALdata <- c(1,1,estInertAL,errInertAL)
	inertAHdata <- c(1,2,estInertAH,errInertAH)
	inertPLdata <- c(2,1,estInertPL,errInertPL)
	inertPHdata <- c(2,2,estInertPH,errInertPH)
		
	plotDataInert <- as.data.frame(rbind(inertALdata, inertAHdata, inertPLdata, inertPHdata))
	names(plotDataInert) <- c("Parameter","Level","sysVar","SE")
	plotDataInert$Parameter <- factor(plotDataInert$Parameter, labels=c("Actor", "Partner"))
	plotDataInert$Level <- factor(plotDataInert$Level, labels=c("Low","High"))
	plotDataInert$errMin <- plotDataInert[ ,3] - plotDataInert[ ,4]
	plotDataInert$errMax <- plotDataInert[ ,3] + plotDataInert[ ,4]
		
	Ymin <- min(basedata$sysVar, na.rm=T)
	Ymax <- max(basedata$sysVar, na.rm=T)
			
	inert <- ggplot(plotDataInert, aes(x=Parameter, y=sysVar, fill=Level)) +
		geom_bar(position=position_dodge(), stat="identity") +
		geom_errorbar(aes(ymin=errMin, ymax=errMax), width=.2, position=position_dodge(.9)) +
		coord_cartesian(ylim=c(Ymin, Ymax)) +
		labs(title="Inertia", subtitle= distName) +
		ylab(sysVarName) +
		theme(plot.title=element_text(size=12))

	coordAData <- basedata
	coordAData$coord0M <- basedata$coord0 - centCoord0[2]
	coordAData$coord1L <- basedata$coord1 - centCoord1[1]
	coordAData$coord0M <- basedata$coord0 - centCoord0[2]
	coordAData$coord1H <- basedata$coord1 - centCoord1[3]
	
	coordPData <- basedata
	coordPData$coord0L <- basedata$coord0 - centCoord0[1]
	coordPData$coord1M <- basedata$coord1 - centCoord1[2]
	coordPData$coord0H <- basedata$coord0 - centCoord0[3]
	coordPData$coord1M <- basedata$coord1 - centCoord1[2]
	
	coordAL <- lm(sysVar ~ coord0M + coord1L, data=coordAData)
	coordAH <- lm(sysVar ~ coord0M + coord1H, data=coordAData)
	coordPL <- lm(sysVar ~ coord0L + coord1M, data=coordPData)
	coordPH <- lm(sysVar ~ coord0H + coord1M, data=coordPData)
	
	estcoordAL <- summary(coordAL)$coefficients[1,1]
	errcoordAL <- summary(coordAL)$coefficients[1,2]
	estcoordAH <- summary(coordAH)$coefficients[1,1]
	errcoordAH <- summary(coordAH)$coefficients[1,2]
	
	estcoordPL <- summary(coordPL)$coefficients[1,1]
	errcoordPL <- summary(coordPL)$coefficients[1,2]
	estcoordPH <- summary(coordPH)$coefficients[1,1]
	errcoordPH <- summary(coordPH)$coefficients[1,2]

	coordALdata <- c(1,1,estcoordAL,errcoordAL)
	coordAHdata <- c(1,2,estcoordAH,errcoordAH)
	coordPLdata <- c(2,1,estcoordPL,errcoordPL)
	coordPHdata <- c(2,2,estcoordPH,errcoordPH)
		
	plotDataCoord <- as.data.frame(rbind(coordALdata, coordAHdata, coordPLdata, coordPHdata))
	names(plotDataCoord) <- c("Parameter","Level","sysVar","SE")
	plotDataCoord$Parameter <- factor(plotDataCoord$Parameter, labels=c("Actor", "Partner"))
	plotDataCoord$Level <- factor(plotDataCoord$Level, labels=c("Low","High"))
	plotDataCoord$errMin <- plotDataCoord[ ,3] - plotDataCoord[ ,4]
	plotDataCoord$errMax <- plotDataCoord[ ,3] + plotDataCoord[ ,4]
		
	Ymin <- min(basedata$sysVar, na.rm=T)
	Ymax <- max(basedata$sysVar, na.rm=T)
			
	coord <- ggplot(plotDataCoord, aes(x=Parameter, y=sysVar, fill=Level)) +
		geom_bar(position=position_dodge(), stat="identity") +
		geom_errorbar(aes(ymin=errMin, ymax=errMax), width=.2, position=position_dodge(.9)) +
		coord_cartesian(ylim=c(Ymin, Ymax)) +
		labs(title="Coordination", subtitle= distName) +
		ylab(sysVarName) +
		theme(plot.title=element_text(size=12))
	
	}
	plotsTemp <- list(inert=inert, coord=coord)
	plots <- gridExtra::grid.arrange(grobs=plotsTemp, nrow=1, ncol=2)

}















