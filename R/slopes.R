

######## This file includes all the functions needed for a patterned-slopes analysis

#' Estimates dyadic linear growth models for each dyad.
#'
#' The observed state variables are predicted from each person's intercept and each person's linear slope over time. Time is centered at the last observation, so intercepts refer to predicted values for each person at the end of the observed time period.
#'
#' @param basedata A dataframe that was produced with the "dataPrep" function.
#' @param dist0name A name for the 0-level of the distinguishing variable (e.g., "Women").
#' @param dist1name A name for the 1-level of the distinguishing variable (e.g., "Men").
#' @param obsName A name for the observed state variables being plotted (e.g., "Emotional Experience").
#' 
#' @return The function returns a list including: 1) the adjusted R^2 for the model for each dyad (called "r2"), 2) the parameter estimates for the model for each dyad (called "paramData", for use in either predicting, or being predicted by, the system variable), and 3) plots of the predicted values against the observed values for each dyad (called "plots"). The plots are also written to the working directory as a pdf file called "linearPlots.pdf"

#' @import ggplot2
#' @export
indivLinear <- function(basedata, dist0name, dist1name, obsName)
{
	newDiD <- unique(factor(basedata$dyad))
	min <- min(basedata$obs, na.rm=T)
	max <- max(basedata$obs, na.rm=T)
	
	r2 <- vector()
	param <- list()
	plots <- list()
	
		for (i in 1:length(newDiD))
		{
			datai <- basedata[basedata$dyad == newDiD[i], ]
			datai$timeEnd <- datai$time - max(datai$time)
			m <- lm(obs ~ -1 + dist0 + dist1 + dist0:timeEnd + dist1:timeEnd , na.action=na.exclude, data=datai)
			r2[[i]] <- summary(m)$adj.r.squared
			param[[i]] <- round(as.numeric(m$coefficients), 5)
			numParam <- length(m$coefficients)
			param[[i]][numParam + 1] <- unique(datai$dyad)
			datai$obsPred <- predict(m)
			datai$role <- factor(datai$dist0, levels=c(0,1), labels=c(dist1name, dist0name)) 
			plotTitle <- as.character(unique(datai$dyad))
						
			plots[[i]] <- ggplot(datai, aes(x=time)) +
				geom_line(aes(y= obs, color=role), linetype="dotted", size= .8, na.rm=T) +
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
		colnames(param) <- c("end0","end1","linear0","linear1","dyad")
		param$endDif <- abs(param$end0 - param$end1)
		temp <- subset(basedata, select=c("id","dyad","sysVar","dist0"))
		temp2 <- unique(temp)
		paramData <- plyr::join(param, temp2)
		
		linearPlots <- gridExtra::marrangeGrob(grobs= plots, ncol=2, nrow=3)
		ggsave('linearPlots.pdf', linearPlots)

	results <- list(r2=r2, paramData=paramData, plots=plots)
}

#' Compares slopes models for predicting the system variable from the dynamic parameters.
#' 
#' The dynamic parameters used in these models come from a set including the linear slope for each person (linear0 and linear1) and the absolute difference between the partner's last observed values (endDif). The 3 models compared are: 1) main effects only (linear0 + linear1 + endDif), 2) the 2-way interaction of slopes, controlling for endDif (linear0*linear1 + endDif), and 3) the 3-way interaction (linear0*linear1*endDif). These are estimated separately for each level of the distinguishing variable. 
#' 
#' @param basedata A dataframe containing the slopes parameter estimates produced by the "indivLinear" function.
#' @param dist0name A name for the 0-level of the distinguishing variable (e.g., "Women").
#' @param dist1name A name for the 1-level of the distinguishing variable (e.g., "Men").
#' @param sysVarName A name for the system variable being predicted (e.g., "Satisfaction").
#' 
#' @return The function returns a list including: 1) the lm objects containing the full results for each model for each level of the distinguishing variable (called "models0" and "models1"), 2) anova output for each model for each level of the distinguishing variable (called "anovas0" and "anovas1"), 3) summary output for each model for each level of the distinguishing variable (called "summaries0" and "summaries1") and 4) adjusted R^2 information for each model for each level of the distinguishing variable (called "adjustR20" and "adjustR21"). The function also displays histograms of the residuals and plots of the predicted values against observed values for each model. 

#' @export
linearSysVarOutCompare <- function(basedata, dist0name, dist1name, sysVarName)
{
	basedata <- basedata[complete.cases(basedata), ] # remove NAs in order to compare nested models	
	basedata0 <- subset(basedata, dist0==1)
	basedata1 <- subset(basedata, dist0==0)

	linear0name <- paste("slope",dist0name, sep="_")
	linear1name <- paste("slope",dist1name, sep="_")
	endDifLinear0name <- paste("endDif",linear0name, sep=":")
	endDifLinear1name <- paste("endDif",linear1name, sep=":")
	linear0linear1name <- paste(linear0name ,linear1name, sep=":")
	endDiflinear0linear1name <- paste("endDif",linear0name ,linear1name, sep=":")

	threeWay0 <- lm(sysVar ~ endDif*linear0*linear1 , data= basedata0)
	names(threeWay0$coefficients) <- c("Intercept","endDif",linear0name,linear1name, endDifLinear0name, endDifLinear1name, linear0linear1name, endDiflinear0linear1name)
	threeWay1 <- lm(sysVar ~ endDif*linear0*linear1 , data= basedata1)
	names(threeWay1$coefficients) <- c("Intercept","endDif",linear0name,linear1name, endDifLinear0name, endDifLinear1name, linear0linear1name, endDiflinear0linear1name)

	twoWay0 <- lm(sysVar ~ endDif + linear0*linear1 , data= basedata0) 
	names(twoWay0$coefficients) <- c("Intercept","endDif",linear0name,linear1name,linear0linear1name)	
	twoWay1 <- lm(sysVar ~ endDif + linear0*linear1 , data= basedata1) 
	names(twoWay1$coefficients) <- c("Intercept","endDif",linear0name,linear1name,linear0linear1name)

	main0 <- lm(sysVar ~ endDif + linear0 + linear1, data= basedata0)
	names(main0$coefficients) <- c("Intercept","endDif",linear0name,linear1name)
	main1 <- lm(sysVar ~ endDif + linear0 + linear1, data= basedata1)
	names(main1$coefficients) <- c("Intercept","endDif",linear0name,linear1name)
	
	main0Anova <- car::Anova(main0, type="III")
	main0Summary <- summary(main0)
	twoWay0Anova <- car::Anova(twoWay0, type="III")
	twoWay0Summary <- summary(twoWay0)
	threeWay0Anova <- car::Anova(threeWay0, type="III")
	threeWay0Summary <- summary(threeWay0)
	
	main1Anova <- car::Anova(main1, type="III")
	main1Summary <- summary(main1)
	twoWay1Anova <- car::Anova(twoWay1, type="III")
	twoWay1Summary <- summary(twoWay1)
	threeWay1Anova <- car::Anova(threeWay1, type="III")
	threeWay1Summary <- summary(threeWay1)
	
	main0Pred <- predict(main0)
	main0R2 <- summary(main0)$adj.r.squared
	
	main1Pred <- predict(main1)
	main1R2 <- summary(main1)$adj.r.squared
	
	twoWay0Pred <- predict(twoWay0)
	twoWay0R2 <- summary(twoWay0)$adj.r.squared
	
	twoWay1Pred <- predict(twoWay1)
	twoWay1R2 <- summary(twoWay1)$adj.r.squared

	threeWay0Pred <- predict(threeWay0)
	threeWay0R2 <- summary(threeWay0)$adj.r.squared
	
	threeWay1Pred <- predict(threeWay1)
	threeWay1R2 <- summary(threeWay1)$adj.r.squared

	min <- min(basedata$sysVar, na.rm=T)
	max <- max(basedata$sysVar, na.rm=T)
	
	ylabName <- paste(sysVarName, "obs", sep="_")
	xlabName <- paste(sysVarName, "pred", sep="_")
	main0name <- paste("MainEffects", dist0name, sep="_")
	main1name <- paste("MainEffects", dist1name, sep="_")
	twoWay0name <- paste("TwoWay", dist0name, sep="_")
	twoWay1name <- paste("TwoWay", dist1name, sep="_")
	threeWay0name <- paste("ThreeWay", dist0name, sep="_")
	threeWay1name <- paste("ThreeWay", dist1name, sep="_")

	par(mfrow=c(2,1))
	hist(residuals(main0))
	plot(basedata0$sysVar ~ main0Pred, xlim=c(min, max), ylim=c(min, max), ylab=ylabName, xlab=xlabName, main=main0name)
	par(mfrow=c(2,1))
	hist(residuals(twoWay0))
	plot(basedata0$sysVar ~ twoWay0Pred, xlim=c(min, max), ylim=c(min, max), ylab=ylabName, xlab=xlabName, main=twoWay0name)
	par(mfrow=c(2,1))
	hist(residuals(threeWay0))
	plot(basedata0$sysVar ~ threeWay0Pred, xlim=c(min, max), ylim=c(min, max), ylab=ylabName, xlab=xlabName, main=threeWay0name)
	
	par(mfrow=c(2,1))
	hist(residuals(main1))
	plot(basedata1$sysVar ~ main1Pred, xlim=c(min, max), ylim=c(min, max), ylab=ylabName, xlab=xlabName, main=main1name)
	par(mfrow=c(2,1))
	hist(residuals(twoWay1))
	plot(basedata1$sysVar ~ twoWay1Pred, xlim=c(min, max), ylim=c(min, max), ylab=ylabName, xlab=xlabName, main=twoWay1name)
	par(mfrow=c(2,1))
	hist(residuals(threeWay1))
	plot(basedata1$sysVar ~ threeWay1Pred, xlim=c(min, max), ylim=c(min, max), ylab=ylabName, xlab=xlabName, main=threeWay1name)
	
	anovas0 <- list(main0Anova=main0Anova, twoWay0Anova=twoWay0Anova, threeWay0Anova=threeWay0Anova)
	anovas1 <- list(main1Anova=main1Anova, twoWay1Anova=twoWay1Anova, threeWay1Anova=threeWay1Anova)
	summaries0 <- list(main0Summary=main0Summary, twoWay0Summary=twoWay0Summary, threeWay0Summary=threeWay0Summary)
	summaries1 <- list(main1Summary=main1Summary, twoWay1Summary=twoWay1Summary, threeWay1Summary=threeWay1Summary)
	r20 <- list(main0R2=main0R2, twoWay0R2=twoWay0R2, threeWay0R2=threeWay0R2)
	r21 <- list(main1R2=main1R2, twoWay1R2=twoWay1R2, threeWay1R2=threeWay1R2)
	models0 <- list(main0=main0, twoWay0=twoWay0, threeWay0=threeWay0)
	models1 <- list(main1=main1, twoWay1=twoWay1, threeWay1=threeWay1)
	
	output <- list(anovas0=anovas0, summaries0=summaries0, adjustR20=r20, anovas1=anovas1, summaries1=summaries1, adjustR21=r21, models0=models0, models1=models1)
}


#' Plots of the system variable predicted from the slopes parameter estimates.
#'
#' Displays bar or line plots for one level of the distinguishing variable, showing model predicted means and standard errors of the system variable for each of the 3 slopes models, which are: 1) at low and high levels of each of the parameter estimates (main effects), 2) at low and high levels of the two slopes, with endDif at an average level (2-way interaction), and 3) at low and high levels of all the parameters (3-way interaction).
#' 
#' @param basedata A dataframe containing the slopes parameter estimates produced by the "indivLinear" function.
#' @param centEnd A vector of low, medium and high centering values for the endDif parameter estimates (e.g., values that indicate typical low, medium and high values of endDif parameter).
#' @param centLine0 A vector of low, medium and high centering values for the slope parameter estimates for the 0-level of the distinguishing variable (e.g., values that indicate typical low, medium and high values of the slope over time for the 0-level partner).
#' @param centLine1 A vector of low, medium and high centering values for the slope parameter estimates for the 1-level of the distinguishing variable (e.g., values that indicate typical low, medium and high values of the slope over time for the 1-level partner).
#' @param sysVarName A name for the system variable being predicted (e.g., "Satisfaction").
#' @param dist A number indicating which level of the distinguishing variable to produce plots for. Must be either 0 or 1.
#' @param dist0name A name for the 0-level of the distinguishing variable (e.g., "Women").
#' @param dist1name A name for the 1-level of the distinguishing variable (e.g., "Men").

#' @import ggplot2
#' @import scales
#' @export

linearSysVarOutPlots <- function(basedata, centEnd, centLine0, centLine1, sysVarName, dist, dist0name, dist1name)
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
	
	mainDataL <- basedata
	mainDataL$endDif <- basedata$endDif - centEnd[1]
	mainDataL$linear0 <- basedata$linear0 - centLine0[1]
	mainDataL$linear1 <- basedata$linear1 - centLine1[1]
	
	mainDataH <- basedata
	mainDataH$endDif <- basedata$endDif - centEnd[3]
	mainDataH$linear0 <- basedata$linear0 - centLine0[3]
	mainDataH$linear1 <- basedata$linear1 - centLine1[3]
	
	endDifL <- lm(sysVar ~ endDif , data= mainDataL, na.action=na.exclude)
	endDifH <- lm(sysVar ~ endDif , data= mainDataH, na.action=na.exclude)
	
	linear0L <- lm(sysVar ~ linear0 , data= mainDataL, na.action=na.exclude)
	linear0H <- lm(sysVar ~ linear0 , data= mainDataH, na.action=na.exclude)
	
	linear1L <- lm(sysVar ~ linear1 , data= mainDataL, na.action=na.exclude)
	linear1H <- lm(sysVar ~ linear1 , data= mainDataH, na.action=na.exclude)
	
	estEndDifL <- summary(endDifL)$coefficients[1,1]
	errEndDifL <- summary(endDifL)$coefficients[1,2]
	estEndDifH <- summary(endDifH)$coefficients[1,1]
	errEndDifH <- summary(endDifH)$coefficients[1,2]
	
	estlinear0L <- summary(linear0L)$coefficients[1,1]
	errlinear0L <- summary(linear0L)$coefficients[1,2]
	estlinear0H <- summary(linear0H)$coefficients[1,1]
	errlinear0H <- summary(linear0H)$coefficients[1,2]
	
	estlinear1L <- summary(linear1L)$coefficients[1,1]
	errlinear1L <- summary(linear1L)$coefficients[1,2]
	estlinear1H <- summary(linear1H)$coefficients[1,1]
	errlinear1H <- summary(linear1H)$coefficients[1,2]
	
	endDifLdata <- c(1, 1, estEndDifL, errEndDifL)
	endDifHdata <- c(1, 2, estEndDifH, errEndDifH)
	linear0Ldata <- c(2, 1, estlinear0L, errlinear0L)
	linear0Hdata <- c(2, 2, estlinear0H, errlinear0H)
	linear1Ldata <- c(3, 1, estlinear1L, errlinear1L)
	linear1Hdata <- c(3, 2, estlinear1H, errlinear1H)
	
	plotData <- as.data.frame(rbind(endDifLdata, endDifHdata, linear0Ldata, linear0Hdata, linear1Ldata, linear1Hdata))
	names(plotData) <- c("Parameter","Level","sysVar","SE")
	slope0name <- paste("slope", dist0name, sep="")
	slope1name <- paste("slope", dist1name, sep="")
	plotData$Parameter <- factor(plotData$Parameter, labels=c("endDif", slope0name, slope1name))
	plotData$Level <- factor(plotData$Level, labels=c("Low","High"))
		
	Ymin <- min(basedata$sysVar, na.rm=T)
	Ymax <- max(basedata$sysVar, na.rm=T)
			
	mainEffects <- ggplot(plotData, aes(x=Parameter, y=sysVar, fill=Level)) +
		geom_bar(position=position_dodge(), stat="identity") +
		geom_errorbar(aes(ymin=sysVar-SE, ymax=sysVar+SE), width=.2, position=position_dodge(.9)) +
		scale_y_continuous(limits=c(Ymin, Ymax), oob=rescale_none) +
		labs(title="MainEffects", subtitle= distName) +
		ylab(sysVarName) +
		theme(plot.title=element_text(size=12))
	
	
	##############
	
	twoWayDataLL <- basedata
	twoWayDataLL$endDif <- basedata$endDif - centEnd[2]
	twoWayDataLL$linear0 <- basedata$linear0 - centLine0[1]
	twoWayDataLL$linear1 <- basedata$linear1 - centLine1[1]
	
	twoWayDataLH <- basedata
	twoWayDataLH$endDif <- basedata$endDif - centEnd[2]
	twoWayDataLH$linear0 <- basedata$linear0 - centLine0[1]
	twoWayDataLH$linear1 <- basedata$linear1 - centLine1[3]
	
	twoWayDataHL <- basedata
	twoWayDataHL$endDif <- basedata$endDif - centEnd[2]
	twoWayDataHL$linear0 <- basedata$linear0 - centLine0[3]
	twoWayDataHL$linear1 <- basedata$linear1 - centLine1[1]
	
	twoWayDataHH <- basedata
	twoWayDataHH$endDif <- basedata$endDif - centEnd[2]
	twoWayDataHH$linear0 <- basedata$linear0 - centLine0[3]
	twoWayDataHH$linear1 <- basedata$linear1 - centLine1[3]
	
	model <- lm(sysVar ~ endDif + linear0*linear1 , data= basedata, na.action=na.exclude)

	LL <- update(model, data=twoWayDataLL)
	estLL <- summary(LL)$coefficients[1,1]
	errLL <- summary(LL)$coefficients[1,2]
	
	LH <- update(model, data=twoWayDataLH)
	estLH <- summary(LH)$coefficients[1,1]
	errLH <- summary(LH)$coefficients[1,2]
	
	HL <- update(model, data=twoWayDataHL)
	estHL <- summary(HL)$coefficients[1,1]
	errHL <- summary(HL)$coefficients[1,2]
	
	HH <- update(model, data=twoWayDataHH)
	estHH <- summary(HH)$coefficients[1,1]
	errHH <- summary(HH)$coefficients[1,2]
	
	LLdata <- c(centLine0[1], centLine1[1], estLL, errLL)
	LHdata <- c(centLine0[1], centLine1[3], estLH, errLH)
	HLdata <- c(centLine0[3], centLine1[1], estHL, errHL)
	HHdata <- c(centLine0[3], centLine1[3], estHH, errHH)
	
	plotData <- as.data.frame(rbind(LLdata, LHdata, HLdata, HHdata))
	names(plotData) <- c(slope0name,slope1name,sysVarName, "SE")
	plotData[ ,1] <- factor(plotData[ ,1])
	plotData[ ,2] <- factor(plotData[ ,2])
	plotData$errMin <- plotData[ ,3] - plotData[ ,4]
	plotData$errMax <- plotData[ ,3] + plotData[ ,4]
	
	Ymin <- min(basedata$sysVar, na.rm=T)
	Ymax <- max(basedata$sysVar, na.rm=T)
				
	twoWay <- ggplot(plotData, aes_string(x=slope0name, y=sysVarName, group=slope1name, color=slope1name)) + 
		geom_line() +
		geom_point() +
		geom_errorbar(aes(ymin=errMin, ymax=errMax), width=.2) + 
		ylim(Ymin, Ymax) +
		ylab(sysVarName) +
		labs(title="TwoWay", subtitle=distName) +
		theme(plot.title=element_text(size=12))
	

	#############
	
	threeWayDataLLL <- basedata
	threeWayDataLLL$endDif <- basedata$endDif - centEnd[1]
	threeWayDataLLL$linear0 <- basedata$linear0 - centLine0[1]
	threeWayDataLLL$linear1 <- basedata$linear1 - centLine1[1]
	
	threeWayDataLLH <- basedata
	threeWayDataLLH$endDif <- basedata$endDif - centEnd[1]
	threeWayDataLLH$linear0 <- basedata$linear0 - centLine0[1]
	threeWayDataLLH$linear1 <- basedata$linear1 - centLine1[3]
	
	threeWayDataLHL <- basedata
	threeWayDataLHL$endDif <- basedata$endDif - centEnd[1]
	threeWayDataLHL$linear0 <- basedata$linear0 - centLine0[3]
	threeWayDataLHL$linear1 <- basedata$linear1 - centLine1[1]
	
	threeWayDataLHH <- basedata
	threeWayDataLHH$endDif <- basedata$endDif - centEnd[1]
	threeWayDataLHH$linear0 <- basedata$linear0 - centLine0[3]
	threeWayDataLHH$linear1 <- basedata$linear1 - centLine1[3]
	
	threeWayDataHLL <- basedata
	threeWayDataHLL$endDif <- basedata$endDif - centEnd[3]
	threeWayDataHLL$linear0 <- basedata$linear0 - centLine0[1]
	threeWayDataHLL$linear1 <- basedata$linear1 - centLine1[1]
	
	threeWayDataHLH <- basedata
	threeWayDataHLH$endDif <- basedata$endDif - centEnd[3]
	threeWayDataHLH$linear0 <- basedata$linear0 - centLine0[1]
	threeWayDataHLH$linear1 <- basedata$linear1 - centLine1[3]
	
	threeWayDataHHL <- basedata
	threeWayDataHHL$endDif <- basedata$endDif - centEnd[3]
	threeWayDataHHL$linear0 <- basedata$linear0 - centLine0[3]
	threeWayDataHHL$linear1 <- basedata$linear1 - centLine1[1]
	
	threeWayDataHHH <- basedata
	threeWayDataHHH$endDif <- basedata$endDif - centEnd[3]
	threeWayDataHHH$linear0 <- basedata$linear0 - centLine0[3]
	threeWayDataHHH$linear1 <- basedata$linear1 - centLine1[3]	
	
	model <- lm(sysVar ~ endDif*linear0*linear1 , data= basedata, na.action=na.exclude)

	LLL <- update(model, data=threeWayDataLLL)
	estLLL <- summary(LLL)$coefficients[1,1]
	errLLL <- summary(LLL)$coefficients[1,2]
	
	LLH <- update(model, data=threeWayDataLLH)
	estLLH <- summary(LLH)$coefficients[1,1]
	errLLH <- summary(LLH)$coefficients[1,2]
	
	LHL <- update(model, data=threeWayDataLHL)
	estLHL <- summary(LHL)$coefficients[1,1]
	errLHL <- summary(LHL)$coefficients[1,2]
	
	LHH <- update(model, data=threeWayDataLHH)
	estLHH <- summary(LHH)$coefficients[1,1]
	errLHH <- summary(LHH)$coefficients[1,2]
	
	HLL <- update(model, data=threeWayDataHLL)
	estHLL <- summary(HLL)$coefficients[1,1]
	errHLL <- summary(HLL)$coefficients[1,2]
	
	HLH <- update(model, data=threeWayDataHLH)
	estHLH <- summary(HLH)$coefficients[1,1]
	errHLH <- summary(HLH)$coefficients[1,2]
	
	HHL <- update(model, data=threeWayDataHHL)
	estHHL <- summary(HHL)$coefficients[1,1]
	errHHL <- summary(HHL)$coefficients[1,2]
	
	HHH <- update(model, data=threeWayDataHHH)
	estHHH <- summary(HHH)$coefficients[1,1]
	errHHH <- summary(HHH)$coefficients[1,2]
	
	LLLdata <- c(centLine0[1], centLine1[1], estLLL, errLLL)
	LLHdata <- c(centLine0[1], centLine1[3], estLLH, errLLH)
	LHLdata <- c(centLine0[3], centLine1[1], estLHL, errLHL)
	LHHdata <- c(centLine0[3], centLine1[3], estLHH, errLHH)
	HLLdata <- c(centLine0[1], centLine1[1], estHLL, errHLL)
	HLHdata <- c(centLine0[1], centLine1[3], estHLH, errHLH)
	HHLdata <- c(centLine0[3], centLine1[1], estHHL, errHHL)
	HHHdata <- c(centLine0[3], centLine1[3], estHHH, errHHH)
		
	plotDataL <- as.data.frame(rbind(LLLdata, LLHdata, LHLdata, LHHdata))
	names(plotDataL) <- c(slope0name,slope1name,sysVarName,"SE")
	plotDataL[ ,1] <- factor(plotDataL[ ,1])
	plotDataL[ ,2] <- factor(plotDataL[ ,2])
	plotDataL$errMin <- plotDataL[ ,3] - plotDataL[ ,4]
	plotDataL$errMax <- plotDataL[ ,3] + plotDataL[ ,4]
	
	plotDataH <- as.data.frame(rbind(HLLdata, HLHdata, HHLdata, HHHdata))
	names(plotDataH) <- c(slope0name,slope1name,sysVarName,"SE")
	plotDataH[ ,1] <- factor(plotDataH[ ,1])
	plotDataH[ ,2] <- factor(plotDataH[ ,2])
	plotDataH$errMin <- plotDataH[ ,3] - plotDataH[ ,4]
	plotDataH$errMax <- plotDataH[ ,3] + plotDataH[ ,4]

	Ymin <- min(basedata$sysVar, na.rm=T)
	Ymax <- max(basedata$sysVar, na.rm=T)
		
	threeWayL <- ggplot(plotDataL, aes_string(slope0name, y=sysVarName, group=slope1name, color=slope1name)) +
		geom_line() +
		geom_point() +
		geom_errorbar(aes(ymin=errMin, ymax=errMax), width=.2) +
		ylim(Ymin, Ymax) +
		ylab(sysVarName) +
		labs(title="ThreeWay, EndDifLow", subtitle= distName) +
		theme(plot.title=element_text(size=12))
	
	threeWayH <- ggplot(plotDataH, aes_string(slope0name, y=sysVarName, group=slope1name, color=slope1name)) +
		geom_line() +
		geom_point() +
		geom_errorbar(aes(ymin=errMin, ymax=errMax), width=.2) +
		ylim(Ymin, Ymax) +
		ylab(sysVarName) +
		labs(title="ThreeWay, EndDifHigh", subtitle= distName) +
		theme(plot.title=element_text(size=12))	
	
	allPlots <- list(mainEffects=mainEffects, twoWay=twoWay, threeWayL=threeWayL, threeWayH=threeWayH)
	plots <- gridExtra::grid.arrange(grobs=allPlots, nrow=2, ncol=2)
}


