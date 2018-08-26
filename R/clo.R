
## This file includes all the functions needed for a coupled oscillator analysis

#### The first two functions were written by Dr. Steven Boker and are available on his website, http://people.virginia.edu/~smb3u/. The method they are implementing is described in the following two publications:

# Boker SM, Nesselroade JR. A method for modeling the intrinsic dynamics of intraindividual variability: Recovering parameters of simulated oscillators in multi-wave panel data. Multivariate Behavioral Research. 2002;37:127–60. 
#  Boker SM, Deboeck PR, Edler C, Keel PK. Generalized local linear approximation of derivatives from time series. In: Chow S, Ferrer E, Hsieh F, editors. Statistical methods for modeling human dynamics. New York: Routledge; 2010. p. 161–78. 

#---------------------------------------------------------
# gllaWMatrix -- Calculates a GLLA linear transformation matrix to 
#                create approximate derivatives
#
# Input:  embed -- Embedding dimension 
#           tau -- Time delay 
#        deltaT -- Interobservation interval
#         order -- Highest order of derivatives (2, 3, or more)

gllaWMatrix <- function(embed, tau, deltaT, order=2) {
    L <- rep(1,embed)
    for(i in 1:order) {
        L <- cbind(L,(((c(1:embed)-mean(1:embed))*tau*deltaT)^i)/factorial(i)) 
    }
    return(L%*%solve(t(L)%*%L))
}

#---------------------------------------------------------
# gllaEmbed -- Creates a time-delay embedding of a variable 
#              given a vector and an optional grouping variable
#              Requires equal interval occasion data ordered by occasion.
#              If multiple individuals, use the ID vector as "groupby"
#
# Input:      x -- vector to embed
#         embed -- Embedding dimension (2 creates an N by 2 embedded matrix)
#           tau -- rows by which to shift x to create each time delay column 
#       groupby -- grouping vector
#         label -- variable label for the columns
#      idColumn -- if TRUE, return ID values in column 1
#                  if FALSE, return the embedding columns only.
#
# Returns:  An embedded matrix where column 1 has the ID values, and the
#           remaining columns are time delay embedded according to the arguments.

gllaEmbed <- function(x, embed, tau, groupby=NA, label="x", idColumn=F) {
    
    minLen <- (tau + 1 + ((embed - 2) * tau))
    if (!is.vector(groupby) | length(groupby[!is.na(groupby[])])<1) {
        groupby <- rep(1,length(x))
    }
    x <- x[!is.na(groupby[])]
    groupby <- groupby[!is.na(groupby[])]
    if (embed < 2 | is.na(embed) | tau < 1 | is.na(tau) | 
        !is.vector(x) | length(x) < minLen)
        return(NA)
    if (length(groupby) != length(x))
        return(NA)
    embeddedMatrix <- matrix(NA, length(x) + (embed*tau), embed+1)
    colNames <- c("ID", paste(label, "0", sep=""))
    for (j in 2:embed) {
        colNames <- c(colNames, paste(label, (j-1)*tau, sep=""))
    }
    dimnames(embeddedMatrix) <- list(NULL, colNames)
    tRow <- 1
    for (i in unique(groupby)) {
        tx <- x[groupby==i]
        if (length(tx) < minLen)
            next
        tLen <- length(tx) - minLen
        embeddedMatrix[tRow:(tRow+tLen), 1] <- i
        for (j in 1:embed) {
            k <- 1 + ((j-1)*tau)
            embeddedMatrix[tRow:(tRow+tLen), j+1] <- tx[k:(k+tLen)]
        }
        tRow <- tRow + tLen + 1
    }
    if (idColumn==TRUE) {
        return(embeddedMatrix[1:(tRow-1),])
    }
    return(embeddedMatrix[1:(tRow-1), 2:(embed+1)])
}

#' Estimates first and second derivatives of an oberved state variable
#'
#' This function makes use of 2 functions written by Steven Boker, "gllaWMatrix" and "gllaEmbed" which are available on his website, http://people.virginia.edu/~smb3u/. It fits an individual linear oscillator model for each person at different combinations of the input parameters (tau, embeds) and returns the input values and period of oscillation that maximize the R^2 for each person. It also estimates first and second derivatives of the observed state variable for each person at the input values that maximize the R^2 for that person and returns a dataframe that contains them.
#'
#' @param basedata A dataframe that was produced with the "dataPrep" function.
#' @param taus A vector containing the values of tau to use. Tau indicates the number of time points to lag in the lagged data matrix (see Boker, S.M., Deboeck, P.R., Edler, C., & Keel, P.K. (2010). Generalized local linear approximation of derivatives from time series. In S.M. Chow & E. Ferrer (Eds.), Statistical Methods for Modeling Human Dynamics: An Interdisciplinary Dialogue (pp. 161-178). New York, NY: Taylor & Francis Group). The first derivative is estimated as the mean of the two adjacent slopes across that number of lags, e.g., if tau = 2 then the estimate of the first derivative at time = t is based on the mean of the slopes left and right of time t across 2 observations each. The second derivative is the difference in the two slopes with respect to time. Tau = 1 is sensitive to noise and increasing its value acts as smoothing. 
#' @param embeds A vector containing the values of embeds to use. Embeds indicates the number of columns in the lagged data matrix. The minimum = 3 for 2nd order derivatives and higher values increase smoothing.
#' @param delta A value indicating the inter-observation interval. For example, if delta = 2, then every second observation is used in the estimation process.
#' 
#' @return The function returns a list including: 1) "data" which is a dataframe containing first and second derivative estimates of an observed state variable, and 2) "fitTable" which shows the maximal R^2 achieved for each person for an individual oscillator model, along with the associated tau, embed and estimated period of oscillation.

#' @export
estDerivs <- function(basedata, taus, embeds, delta)
{
   basedata <- basedata[complete.cases(basedata), ] 
   params <- expand.grid(taus=taus, embeds=embeds)
   newId <- unique(factor(basedata$id))

   derivData <- list()
   fitTable <- list()

	for(p in 1:length(newId))
	{
	   r <- list()
	   datai <- basedata[basedata$id == newId[p],]

		for(i in 1:nrow(params))
		{
		   # Estimate derivatives for each parameter combination
		   obsEmbed <- gllaEmbed(datai$obs_deTrend, tau=params[i,1], embed=params[i,2])
		   p_obsEmbed <- gllaEmbed(datai$p_obs_deTrend, tau=params[i,1], embed=params[i,2])   

		   obsLLA <- gllaWMatrix(tau=params[i,1], embed=params[i,2], deltaT=delta, order=2)
		   obsDeriv <- obsEmbed[,1:dim(obsEmbed)[2]] %*% obsLLA
		   p_obsDeriv <- p_obsEmbed[,1:dim(p_obsEmbed)[2]] %*% obsLLA

		   # combine actor/partner derivatives and name columns
		   deriv <- cbind(obsDeriv, p_obsDeriv)
		   dimnames(deriv) <- list(NULL, c("obs_deTrend","d1","d2","p_obs_deTrend","p_d1","p_d2"))
		   deriv <- as.data.frame(deriv)

		   # fit CLO and get R^2 for that combination of parameters
		   out <- lm(d2 ~ obs_deTrend + d1 + p_obs_deTrend + p_d1, data=deriv)
		   r[i] <- summary(out)$adj.r.squared 
		}
		
	   ## get maximum R^2 for given person and the corresponding tau & embed values
	   maxR <- max(unlist(r))
	   paramRow <- which(r==maxR)
	   select <- params[paramRow,]
	   tau <- select[1,1]
	   embed <- select[1,2]

	   ## get period (1 cycle (peak to peak)) in given time units associated with highest R^2

	   if (out$coefficients[2] >= 0) {print("error: frequency parameter is not negative")}
	   n <- abs(as.numeric(out$coefficients[2]))
	   period <- (2*pi) / (sqrt(n))

	   # Estimate derivatives with selected tau and embed for each person
	   obsEmbed <- gllaEmbed(datai$obs_deTrend, tau=tau, embed=embed)
	   p_obsEmbed <- gllaEmbed(datai$p_obs_deTrend, tau=tau, embed=embed)   

	   obsLLA <- gllaWMatrix(tau=tau, embed=embed, deltaT=delta, order=2)
	   obsDeriv <- obsEmbed[,1:dim(obsEmbed)[2]] %*% obsLLA
	   p_obsDeriv <- p_obsEmbed[,1:dim(p_obsEmbed)[2]] %*% obsLLA

	   # add time, recover original ids, distinguisher and moderator variables, add to derivatives and 
	   # name columns
	   idLength <- dim(obsDeriv)[1]
	   time <- seq_len(idLength)
	   dyad <- rep(unique(datai$dyad, idLength))
	   id <- rep(unique(datai$id), idLength)
	   dist0 <- rep(unique(datai$dist0), idLength)
	   dist1 <- rep(unique(datai$dist1), idLength)
	   sysVar <- rep(unique(datai$sysVar), idLength)
	   deriv <- cbind(dyad, id, time, obsDeriv, p_obsDeriv, dist0, dist1, sysVar)
	   dimnames(deriv) <- list(NULL,c("dyad","id","time","obs_deTrend","d1","d2","p_obs_deTrend","p_d1","p_d2","dist0", 
	   "dist1", "sysVar"))
	   deriv <- as.data.frame(deriv)

	   derivData[[p]] <- deriv
	   fitTable[[p]] <- c("id"= unique(deriv$id), "tau"=tau, "embed"= embed, 
	   "Rsqr"= maxR, "Period"=period)
	}

   ## output, which includes the derivative data and the fit table

   derivD <- as.data.frame(do.call(rbind, derivData))
   fitTable <- matrix(unlist(fitTable), nrow= length(newId), byrow=T)
   dimnames(fitTable) <- list(NULL, c("id","tau","embed","Max-Rsqr","Period"))
   fitTable <- round(fitTable, 2)

   derivOut <- list("data"=derivD, "fitTable"= fitTable)
}

#' Estimates a "coupled oscillator" model for each dyad.
#' 
#' The second derivatives of the observed state variables (with linear trends removed) are predicted from each person's own and partner's observed state variables (again with linear trends removed), as well as each person's own and partner's first derivatives of the observed state variables (again with linear trends removed).
#'
#' @param basedata A dataframe that was produced with the "estDerivs" function.
#' @param dist0name A name for the 0-level of the distinguishing variable (e.g., "Women").
#' @param dist1name A name for the 1-level of the distinguishing variable (e.g., "Men").
#' @param obsName A name for the observed state variables being plotted (e.g., "Emotional Experience").
#' 
#' @return The function returns a list including: 1) the adjusted R^2 for the model for each dyad (called "r2"), 2) the parameter estimates for the model for each dyad (called "paramData", for use in either predicting, or being predicted by, the system variable), and 3) plots of the predicted values against the observed values for each dyad (called "plots"). The plots are also written to the working directory as a pdf file called "cloPlotsCouple.pdf"

#' @import ggplot2
#' @export
indivCloCouple <- function(basedata, dist0name, dist1name, obsName)
{
	newDiD <- unique(factor(basedata$dyad))
	basedata <- basedata[complete.cases(basedata), ]
	min <- min(basedata$obs_deTrend, na.rm=T)
	max <- max(basedata$obs_deTrend, na.rm=T)
	
	r2 <- vector()
	param <- list()
	plots <- list()
	
		for (i in 1:length(newDiD))
		{
			statedatai <- basedata[basedata$dyad == newDiD[i] & basedata$dist0 == 1,] 
  			maxtime <- max(statedatai$time) 
 			plotTimes <- seq(1, maxtime, by=1)
 			start <- suppressWarnings(subset(statedatai, time==c(1:5), select=c(obs_deTrend, p_obs_deTrend)))
			y1 <- mean(start$obs_deTrend, na.rm=T)
 			y2 <- 0
			y3 <- mean(start$p_obs_deTrend, na.rm=T)
 			y4 <- 0
 			statei <- c("y1"=y1, "y2"=y2, "y3"=y3, "y4"=y4)
			
			datai <- basedata[basedata$dyad == newDiD[i], ]
			m <- lm(d2 ~ dist0:obs_deTrend + dist0:d1 + dist0:p_obs_deTrend + dist0:p_d1 + dist1:obs_deTrend + dist1:d1 + dist1:p_obs_deTrend + dist1:p_d1 -1, na.action=na.exclude, data=datai)
			r2[[i]] <- summary(m)$adj.r.squared
			param[[i]] <- round(as.numeric(m$coefficients), 5)
			numParam <- length(m$coefficients)
			param[[i]][numParam + 1] <- unique(datai$dyad)
			
			obs_0 <- param[[i]][1]	
			d1_0 <- param[[i]][2]
			p_obs_0 <- param[[i]][3]
			p_d1_0 <- param[[i]][4]
			obs_1 <- param[[i]][5]
			d1_1 <- param[[i]][6]
			p_obs_1 <- param[[i]][7]
			p_d1_1 <- param[[i]][8]
			paramClo <- c("obs_0"= obs_0, "d1_0"= d1_0, "p_obs_0"= p_obs_0, "p_d1_0"=p_d1_0, "obs_1"=obs_1, "d1_1"=d1_1, "p_obs_1"= p_obs_1, "p_d1_1"= p_d1_1)
			
			temp <- as.data.frame(deSolve::ode(y=statei, times=plotTimes, func=cloCoupleOde, parms= paramClo))
			temp2 <- subset(temp, select=-c(y2, y4))
			names(temp2) <- c("time","d0.pred","d1.pred")
			temp2$dyad <- statedatai$dyad
			temp3 <- reshape(temp2, direction='long', varying=c("d0.pred","d1.pred"), timevar="role", times=c("d0","d1"), v.names=c("pred"), idvar="time")
			temp3$id <- ifelse(temp3$role == "d0", temp3$dyad, temp3$dyad+500)
			temp4 <- suppressMessages(plyr::join(datai, temp3))
			temp4$roleNew <- factor(temp4$role, levels=c("d0","d1"), labels=c(dist0name, dist1name)) 
			
			plotData <- temp4[complete.cases(temp4), ]	
			plotTitle <- as.character(unique(datai$dyad))
						
			plots[[i]] <- ggplot(plotData, aes(x=time)) +
				geom_line(aes(y= obs_deTrend, color=roleNew), linetype="dotted", size= .8, na.rm=T) +
				geom_line(aes(y=pred, color=roleNew), size= .8, na.rm=T) + 
				scale_color_manual(name="Role", values=c("red","blue")) +
				ylab(obsName) +
				ylim(min, max) +
				annotate("text", x=-Inf, y=-Inf, hjust=0, vjust=0, label="Dots = Observed; Lines = Predicted", size=3) +
				labs(title= "Dyad ID:", subtitle= plotTitle) +
				theme(plot.title=element_text(size=11)) +
				theme(plot.subtitle=element_text(size=10))			
		}
	
		param <- as.data.frame(do.call(rbind, param))
		colnames(param) <- c("obs_0","d1_0","p_obs_0","p_d1_0","obs_1","d1_1","p_obs_1","p_d1_1","dyad")
		temp <- subset(basedata, select=c("id","dyad","sysVar","dist0"))
		temp2 <- unique(temp)
		paramData <- suppressMessages(plyr::join(param, temp2))
		
		cloPlots <- gridExtra::marrangeGrob(grobs= plots, ncol=2, nrow=3)
		ggsave('cloPlotsCouple.pdf', cloPlots)

	results <- list(r2=r2, paramData=paramData, plots=plots)
}


#############################

#' Estimates an "un-coupled oscillator" model (e.g., only self frequency and damping with no coupling) for each dyad.
#' 
#' The second derivatives of the observed state variables (with linear trends removed) are predicted from each person's own observed state variable (again with linear trends removed), as well as each person's own first derivatives of the observed state variables (again with linear trends removed).
#'
#' @param basedata A dataframe that was produced with the "estDerivs" function.
#' @param dist0name A name for the 0-level of the distinguishing variable (e.g., "Women").
#' @param dist1name A name for the 1-level of the distinguishing variable (e.g., "Men").
#' @param obsName A name for the observed state variables being plotted (e.g., "Emotional Experience").
#' 
#' @return The function returns a list including: 1) the adjusted R^2 for the model for each dyad (called "r2"), 2) the parameter estimates for the model for each dyad (called "paramData", for use in either predicting, or being predicted by, the system variable), and 3) plots of the predicted values against the observed values for each dyad (called "plots"). The plots are also written to the working directory as a pdf file called "cloPlotsUncouple.pdf"

#' @import ggplot2
#' @export
indivCloUncouple <- function(basedata, dist0name, dist1name, obsName)
{
	newDiD <- unique(factor(basedata$dyad))
	basedata <- basedata[complete.cases(basedata), ]
	min <- min(basedata$obs_deTrend, na.rm=T)
	max <- max(basedata$obs_deTrend, na.rm=T)
	
	r2 <- vector()
	param <- list()
	plots <- list()
	
		for (i in 1:length(newDiD))
		{
			statedatai <- basedata[basedata$dyad == newDiD[i] & basedata$dist0 == 1,] 
  			maxtime <- max(statedatai$time) 
 			plotTimes <- seq(1, maxtime, by=1)
 			start <- suppressWarnings(subset(statedatai, time==c(1:5), select=c(obs_deTrend, p_obs_deTrend)))
			y1 <- mean(start$obs_deTrend, na.rm=T)
 			y2 <- 0
			y3 <- mean(start$p_obs_deTrend, na.rm=T)
 			y4 <- 0
 			statei <- c("y1"=y1, "y2"=y2, "y3"=y3, "y4"=y4)
			
			datai <- basedata[basedata$dyad == newDiD[i], ]
			m <- lm(d2 ~ dist0:obs_deTrend + dist0:d1 + dist1:obs_deTrend + dist1:d1 -1, na.action=na.exclude, data=datai)
			r2[[i]] <- summary(m)$adj.r.squared
			param[[i]] <- round(as.numeric(m$coefficients), 5)
			numParam <- length(m$coefficients)
			param[[i]][numParam + 1] <- unique(datai$dyad)
			
			obs_0 <- param[[i]][1]	
			d1_0 <- param[[i]][2]
			obs_1 <- param[[i]][3]
			d1_1 <- param[[i]][4]
			paramClo <- c("obs_0"= obs_0, "d1_0"= d1_0, "obs_1"=obs_1, "d1_1"=d1_1)
			
			temp <- as.data.frame(deSolve::ode(y=statei, times=plotTimes, func=cloUncoupleOde, parms= paramClo))
			temp2 <- subset(temp, select=-c(y2, y4))
			names(temp2) <- c("time","d0.pred","d1.pred")
			temp2$dyad <- statedatai$dyad
			temp3 <- reshape(temp2, direction='long', varying=c("d0.pred","d1.pred"), timevar="role", times=c("d0","d1"), v.names=c("pred"), idvar="time")
			temp3$id <- ifelse(temp3$role == "d0", temp3$dyad, temp3$dyad+500)
			temp4 <- suppressMessages(plyr::join(datai, temp3))
			temp4$roleNew <- factor(temp4$role, levels=c("d0","d1"), labels=c(dist0name, dist1name)) 
			
			plotData <- temp4[complete.cases(temp4), ]	
			plotTitle <- as.character(unique(datai$dyad))
						
			plots[[i]] <- ggplot(plotData, aes(x=time)) +
				geom_line(aes(y= obs_deTrend, color=roleNew), linetype="dotted", size= .8, na.rm=T) +
				geom_line(aes(y=pred, color=roleNew), size= .8, na.rm=T) + 
				scale_color_manual(name="Role", values=c("red","blue")) +
				ylab(obsName) +
				ylim(min, max) +
				annotate("text", x=-Inf, y=-Inf, hjust=0, vjust=0, label="Dots = Observed; Lines = Predicted", size=3) +
				labs(title= "Dyad ID:", subtitle= plotTitle) +
				theme(plot.title=element_text(size=11)) +
				theme(plot.subtitle=element_text(size=10))			
		}
	
		param <- as.data.frame(do.call(rbind, param))
		colnames(param) <- c("obs_0","d1_0","obs_1","d1_1","dyad")
		temp <- subset(basedata, select=c("id","dyad","sysVar","dist0"))
		temp2 <- unique(temp)
		paramData <- suppressMessages(plyr::join(param, temp2))
		
		cloPlots <- gridExtra::marrangeGrob(grobs= plots, ncol=2, nrow=3)
		ggsave('cloPlotsUncouple.pdf', cloPlots)

	results <- list(r2=r2, paramData=paramData, plots=plots)
}


#' Provides the equation for a coupled oscillator model for the differential equation solver (ode) to plot

cloCoupleOde <- function(t, state, parameters)
{
	with(as.list(c(state, parameters)), 
	{
		dy1 <- y2
		dy2 <- y1*obs_0 + y2*d1_0 + y3*p_obs_0 + y4*p_d1_0
		dy3 <- y4
		dy4 <- y3*obs_1 + y4*d1_1 + y1*p_obs_1 + y2*p_d1_1
		list(c(dy1, dy2, dy3, dy4))		
	})
}


#' Provides the equation for an un-coupled oscillator model for the differential equation solver (ode) to plot

cloUncoupleOde <- function(t, state, parameters)
{
	with(as.list(c(state, parameters)), 
	{
		dy1 <- y2
		dy2 <- y1*obs_0 + y2*d1_0 
		dy3 <- y4
		dy4 <- y3*obs_1 + y4*d1_1 
		list(c(dy1, dy2, dy3, dy4))		
	})
}


#' Compares model fit for the uncoupled and coupled oscillator for each dyad's state trajectories using an R-square comparison. 
#' 
#' Fits an uncoupled and coupled oscillator model to each dyad's observed state variables and returns the adjusted R-squares, along with the difference between them (coupled - uncoupled, so positive values indicate better fit for the more complex model).
#'
#' @param basedata A dataframe that was produced with the "estDerivs" function.
#' 
#' @return The function returns a named list including: 1) the adjusted R^2 for the uncoupled model for each dyad (called "r2uncouple"), 2) the adjusted R^2 for the coupled model for each dyad (called "r2couple"), and 3) the difference between the R-squares for each dyad (coupled - uncoupled, called "r2dif").

#' @export
cloIndivCompare <- function(basedata)
{
	newDiD <- unique(factor(basedata$dyad))
	basedata <- basedata[complete.cases(basedata), ]

	r2uncouple <- vector()
	r2couple <- vector()
	r2dif <- vector()
	
		for (i in 1:length(newDiD))
		{
			datai <- basedata[basedata$dyad == newDiD[i], ]
			uncouple <- lm(d2 ~ dist0:obs_deTrend + dist0:d1 + dist1:obs_deTrend + dist1:d1 -1, na.action=na.exclude, data=datai)
			couple <- lm(d2 ~ dist0:obs_deTrend + dist0:d1 + dist0:p_obs_deTrend + dist0:p_d1 + dist1:obs_deTrend + dist1:d1 + dist1:p_obs_deTrend + dist1:p_d1-1, na.action=na.exclude, data=datai)	
			r2uncouple[[i]] <- summary(uncouple)$adj.r.squared
			r2couple[[i]] <- summary(couple)$adj.r.squared	
			r2dif[[i]] <- r2couple[[i]] - r2uncouple[[i]]
		}			
		output <- list(r2uncouple=r2uncouple, r2couple=r2couple, r2dif=r2dif)
		}


#' Compares coupled and uncoupled clo models for predicting the system variable from the dynamic parameters.
#' 
#' The dynamic parameters used in these models come from a set including people's own and partner's influence on their frequency (freqS0, freqS1, freqP0, freqP1) and own and partner's influence on their damping/amplification (dampS0, dampS1, dampP0, dampP1). The 2 models compared are the uncoupled, self-only CLO (freqS + dampS) and the full coupled CLO (freqS + dampS + freqP + dampP). These are estimated separately for each level of the distinguishing variable. 
#' 
#' @param basedata A dataframe containing the coupled oscillator parameter estimates produced by the "indivClo" function.
#' @param dist0name A name for the 0-level of the distinguishing variable (e.g., "Women").
#' @param dist1name A name for the 1-level of the distinguishing variable (e.g., "Men").
#' @param sysVarName A name for the system variable being predicted (e.g., "Satisfaction").
#' 
#' @return The function returns a list including: 1) the lm objects containing the full results for each model for each level of the distinguishing variable (called "models0" and "models1"), 2) anova output for each model for each level of the distinguishing variable (called "anovas0" and "anovas1"), 3) summary output for each model for each level of the distinguishing variable (called "summaries0" and "summaries1") and 4) adjusted R^2 information for each model for each level of the distinguishing variable (called "adjustR20" and "adjustR21"). The function also displays histograms of the residuals and plots of the predicted values against observed values for each model. 

#' @export
cloSysVarOutCompare <- function(basedata, dist0name, dist1name, sysVarName)
{
	basedata <- basedata[complete.cases(basedata), ] 
	basedata0 <- subset(basedata, dist0==1)
	basedata1 <- subset(basedata, dist0==0)

	freqSname <- paste("freq","self", sep="_")
	dampSname <- paste("damp","self", sep="_")
	freqPname <- paste("freq","partner", sep="_")
	dampPname <- paste("damp","partner", sep="_")

	uncoupled0 <- lm(sysVar ~ obs_0 + d1_0 , data= basedata0)
	names(uncoupled0$coefficients) <- c("intercept", freqSname, dampSname)
	uncoupled1 <- lm(sysVar ~ obs_1 + d1_1 , data= basedata1)
	names(uncoupled1$coefficients) <- c("intercept", freqSname, dampSname)
	
	couple0 <- lm(sysVar ~ obs_0 + d1_0 + p_obs_0 + p_d1_0, data= basedata0)
	names(couple0$coefficients) <- c("intercept", freqSname, dampSname, freqPname, dampPname)
	couple1 <- lm(sysVar ~ obs_1 + d1_1 + p_obs_1 + p_d1_1, data= basedata1)
	names(couple1$coefficients) <- c("intercept", freqSname, dampSname, freqPname, dampPname)
	
	uncoupled0Anova <- car::Anova(uncoupled0, type="III")
	uncoupled0Summary <- summary(uncoupled0)
	couple0Anova <- car::Anova(couple0, type="III")
	couple0Summary <- summary(couple0)
	
	uncoupled1Anova <- car::Anova(uncoupled1, type="III")
	uncoupled1Summary <- summary(uncoupled1)
	couple1Anova <- car::Anova(couple1, type="III")
	couple1Summary <- summary(couple1)

	uncoupled0Pred <- predict(uncoupled0)
	uncoupled0R2 <- summary(uncoupled0)$adj.r.squared
	couple0Pred <- predict(couple0)
	couple0R2 <- summary(couple0)$adj.r.squared
	
	uncoupled1Pred <- predict(uncoupled1)
	uncoupled1R2 <- summary(uncoupled1)$adj.r.squared
	couple1Pred <- predict(couple1)
	couple1R2 <- summary(couple1)$adj.r.squared

	min0 <- min(basedata0$sysVar, na.rm=T)
	max0 <- max(basedata0$sysVar, na.rm=T)
	min1 <- min(basedata1$sysVar, na.rm=T)
	max1 <- max(basedata1$sysVar, na.rm=T)

	ylabName <- paste(sysVarName, "obs", sep=" ")
	xlabName <- paste(sysVarName, "pred", sep=" ")
	uncoupled0name <- paste("Uncoupled (Self Only)", dist0name, sep=" ")
	uncoupled1name <- paste("Uncoupled (Self Only)", dist1name, sep="_")
	couple0name <- paste("Coupled (Self & Partner)", dist0name, sep=" ")
	couple1name <- paste("Coupled (Self & Partner)", dist1name, sep=" ")

	par(mfrow=c(2,1))
	hist(residuals(uncoupled0))
	plot(basedata0$sysVar ~ uncoupled0Pred, xlim=c(min0, max0), ylim=c(min0, max0), ylab=ylabName, xlab=xlabName, main=uncoupled0name)
	par(mfrow=c(2,1))
	hist(residuals(couple0))
	plot(basedata0$sysVar ~ couple0Pred, xlim=c(min0, max0), ylim=c(min0, max0), ylab=ylabName, xlab=xlabName, main=couple0name)
	
	par(mfrow=c(2,1))
	hist(residuals(uncoupled1))
	plot(basedata1$sysVar ~ uncoupled1Pred, xlim=c(min0, max0), ylim=c(min0, max0), ylab=ylabName, xlab=xlabName, main=uncoupled1name)
	par(mfrow=c(2,1))
	hist(residuals(couple1))
	plot(basedata1$sysVar ~ couple1Pred, xlim=c(min0, max0), ylim=c(min0, max0), ylab=ylabName, xlab=xlabName, main=couple1name)

	models0 <- list(uncoupled0=uncoupled0, couple0=couple0)
	models1 <- list(uncoupled1=uncoupled1, couple1=couple1)
	anovas0 <- list(uncoupled0Anova=uncoupled0Anova, couple0Anova=couple0Anova)
	anovas1 <- list(uncoupled1Anova=uncoupled1Anova, couple1Anova=couple1Anova)
	summaries0 <- list(uncoupled0Summary=uncoupled0Summary, couple0Summary=couple0Summary)
	summaries1 <- list(uncoupled1Summary=uncoupled1Summary, couple1Summary=couple1Summary)	
	r20 <- list(uncoupled0R2=uncoupled0R2, couple0R2=couple0R2)
	r21 <- list(uncoupled1R2=uncoupled1R2, couple1R2=couple1R2)
	
	output <- list(models0=models0, models1=models1, anovas0=anovas0, summaries0=summaries0, adjustR20=r20, anovas1=anovas1, summaries1=summaries1, adjustR21=r21)
}

#' Plots of the system variable predicted from the uncoupled version of the clo model.
#'
#' Displays bar plots for one level of the distinguishing variable, showing model predicted means and standard errors of the system variable at low and high levels of a person's own frequency and damping/amplification parameter estimates.
#' 
#' @param basedata A dataframe containing the clo parameter estimates produced by the "indivClo" function.
#' @param centFreqSelf0 A vector of low, medium and high centering values for the frequency parameter estimates for the 0-level of the distinguishing variable (e.g., values that indicate typical low, medium and high values of frequency for the 0-level partner).
#' @param centDampSelf0 A vector of low, medium and high centering values for the damping/amplification parameter estimates for the 0-level of the distinguishing variable (e.g., values that indicate typical low, medium and high values of damping/amplification for the 0-level partner).
#' @param centFreqPart0 A vector of low, medium and high centering values for the partner influence on frequency parameter estimates for the 0-level of the distinguishing variable (e.g., values that indicate typical low, medium and high values of partner influence on frequency for the 0-level partner).
#' @param centDampPart0 A vector of low, medium and high centering values for the partner influence on damping/amplification parameter estimates for the 0-level of the distinguishing variable (e.g., values that indicate typical low, medium and high values of partner influence on damping/amplification for the 0-level partner).
#' @param centFreqSelf1 A vector of low, medium and high centering values for the frequency parameter estimates for the 1-level of the distinguishing variable (e.g., values that indicate typical low, medium and high values of frequency for the 1-level partner).
#' @param centDampSelf1 A vector of low, medium and high centering values for the damping/amplification parameter estimates for the 1-level of the distinguishing variable (e.g., values that indicate typical low, medium and high values of damping/amplification for the 1-level partner).
#' @param centFreqPart1 A vector of low, medium and high centering values for the partner influence on frequency parameter estimates for the 1-level of the distinguishing variable (e.g., values that indicate typical low, medium and high values of partner influence on frequency for the 1-level partner).
#' @param centDampPart1 A vector of low, medium and high centering values for the partner influence on damping/amplification parameter estimates for the 1-level of the distinguishing variable (e.g., values that indicate typical low, medium and high values of partner influence on damping/amplification for the 1-level partner).
#' @param sysVarName A name for the system variable being predicted (e.g., "Satisfaction").
#' @param dist A number indicating which level of the distinguishing variable to produce plots for. Must be either 0 or 1.
#' @param dist0name A name for the 0-level of the distinguishing variable (e.g., "Women").
#' @param dist1name A name for the 1-level of the distinguishing variable (e.g., "Men").

#' @import ggplot2
#' @export
uncoupleSysVarOutPlots <- function(basedata, centFreqSelf0, centDampSelf0, centFreqPart0, centDampPart0, centFreqSelf1, centDampSelf1, centFreqPart1, centDampPart1, sysVarName, dist, dist0name, dist1name)
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
	freqData <- basedata
	freqData$freqL <- basedata$obs_0 - centFreqSelf0[1]
	freqData$freqH <- basedata$obs_0 - centFreqSelf0[3]
	freqData$dampM <- basedata$d1_0 - centDampSelf0[2]
	
	dampData <- basedata
	dampData$dampL <- basedata$d1_0 - centDampSelf0[1]
	dampData$dampH <- basedata$d1_0 - centDampSelf0[3]
	dampData$freqM <- basedata$obs_0 - centFreqSelf0[2]
	
	freqL <- lm(sysVar ~ freqL + dampM, data=freqData)
	freqH <- lm(sysVar ~ freqH + dampM, data=freqData)
	dampL <- lm(sysVar ~ freqM + dampL, data=dampData)
	dampH <- lm(sysVar ~ freqM + dampH, data=dampData)
	
	estFreqL <- summary(freqL)$coefficients[1,1]
	errFreqL <- summary(freqL)$coefficients[1,2]
	estFreqH <- summary(freqH)$coefficients[1,1]
	errFreqH <- summary(freqH)$coefficients[1,2]
	
	estDampL <- summary(dampL)$coefficients[1,1]
	errDampL <- summary(dampL)$coefficients[1,2]
	estDampH <- summary(dampH)$coefficients[1,1]
	errDampH <- summary(dampH)$coefficients[1,2]

	freqLdata <- c(1,1,estFreqL,errFreqL)
	freqHdata <- c(1,2,estFreqH,errFreqH)
	dampLdata <- c(2,1,estDampL,errDampL)
	dampHdata <- c(2,2,estDampH,errDampH)
		
	plotData <- as.data.frame(rbind(freqLdata, freqHdata, dampLdata, dampHdata))
	names(plotData) <- c("Parameter","Level","sysVar","SE")
	plotData$Parameter <- factor(plotData$Parameter, labels=c("Freq", "Damp"))
	plotData$Level <- factor(plotData$Level, labels=c("Low","High"))
	plotData$errMin <- plotData[ ,3] - plotData[ ,4]
	plotData$errMax <- plotData[ ,3] + plotData[ ,4]
		
	Ymin <- min(basedata$sysVar, na.rm=T)
	Ymax <- max(basedata$sysVar, na.rm=T)
			
	self <- ggplot(plotData, aes(x=Parameter, y=sysVar, fill=Level)) +
		geom_bar(position=position_dodge(), stat="identity") +
		geom_errorbar(aes(ymin=errMin, ymax=errMax), width=.2, position=position_dodge(.9)) +
		coord_cartesian(ylim=c(Ymin, Ymax)) +
		labs(title="Uncoupled (Self-Only)", subtitle= distName) +
		ylab(sysVarName) +
		theme(plot.title=element_text(size=12))
	
	} else if (dist==1)
	
	{
	freqData <- basedata
	freqData$freqL <- basedata$obs_1 - centFreqSelf1[1]
	freqData$freqH <- basedata$obs_1 - centFreqSelf1[3]
	freqData$dampM <- basedata$d1_1 - centDampSelf1[2]
	
	dampData <- basedata
	dampData$dampL <- basedata$d1_1 - centDampSelf1[1]
	dampData$dampH <- basedata$d1_1 - centDampSelf1[3]
	dampData$freqM <- basedata$obs_1 - centFreqSelf1[2]
	
	freqL <- lm(sysVar ~ freqL + dampM, data=freqData)
	freqH <- lm(sysVar ~ freqH + dampM, data=freqData)
	dampL <- lm(sysVar ~ freqM + dampL, data=dampData)
	dampH <- lm(sysVar ~ freqM + dampH, data=dampData)
	
	estFreqL <- summary(freqL)$coefficients[1,1]
	errFreqL <- summary(freqL)$coefficients[1,2]
	estFreqH <- summary(freqH)$coefficients[1,1]
	errFreqH <- summary(freqH)$coefficients[1,2]
	
	estDampL <- summary(dampL)$coefficients[1,1]
	errDampL <- summary(dampL)$coefficients[1,2]
	estDampH <- summary(dampH)$coefficients[1,1]
	errDampH <- summary(dampH)$coefficients[1,2]

	freqLdata <- c(1,1,estFreqL,errFreqL)
	freqHdata <- c(1,2,estFreqH,errFreqH)
	dampLdata <- c(2,1,estDampL,errDampL)
	dampHdata <- c(2,2,estDampH,errDampH)
		
	plotData <- as.data.frame(rbind(freqLdata, freqHdata, dampLdata, dampHdata))
	names(plotData) <- c("Parameter","Level","sysVar","SE")
	plotData$Parameter <- factor(plotData$Parameter, labels=c("Freq", "Damp"))
	plotData$Level <- factor(plotData$Level, labels=c("Low","High"))
	plotData$errMin <- plotData[ ,3] - plotData[ ,4]
	plotData$errMax <- plotData[ ,3] + plotData[ ,4]
		
	Ymin <- min(basedata$sysVar, na.rm=T)
	Ymax <- max(basedata$sysVar, na.rm=T)
			
	self <- ggplot(plotData, aes(x=Parameter, y=sysVar, fill=Level)) +
		geom_bar(position=position_dodge(), stat="identity") +
		geom_errorbar(aes(ymin=errMin, ymax=errMax), width=.2, position=position_dodge(.9)) +
		coord_cartesian(ylim=c(Ymin, Ymax)) +
		labs(title="Uncoupled (Self-Only)", subtitle= distName) +
		ylab(sysVarName) +
		theme(plot.title=element_text(size=12))
	
	}
	plots <- self
}

#' Plots of the system variable predicted from the full coupled version of the clo model.
#'
#' Displays bar plots for one level of the distinguishing variable, showing model predicted means and standard errors of the system variable at low and high levels of a person's own frequency and damping/amplification parameter estimates and low and high levels of the partner's influence on their frequency and damping/amplification parameter estimates.
#' 
#' @param basedata A dataframe containing the clo parameter estimates produced by the "indivClo" function.
#' @param centFreqSelf0 A vector of low, medium and high centering values for the frequency parameter estimates for the 0-level of the distinguishing variable (e.g., values that indicate typical low, medium and high values of frequency for the 0-level partner).
#' @param centDampSelf0 A vector of low, medium and high centering values for the damping/amplification parameter estimates for the 0-level of the distinguishing variable (e.g., values that indicate typical low, medium and high values of damping/amplification for the 0-level partner).
#' @param centFreqPart0 A vector of low, medium and high centering values for the partner influence on frequency parameter estimates for the 0-level of the distinguishing variable (e.g., values that indicate typical low, medium and high values of partner influence on frequency for the 0-level partner).
#' @param centDampPart0 A vector of low, medium and high centering values for the partner influence on damping/amplification parameter estimates for the 0-level of the distinguishing variable (e.g., values that indicate typical low, medium and high values of partner influence on damping/amplification for the 0-level partner).
#' @param centFreqSelf1 A vector of low, medium and high centering values for the frequency parameter estimates for the 1-level of the distinguishing variable (e.g., values that indicate typical low, medium and high values of frequency for the 1-level partner).
#' @param centDampSelf1 A vector of low, medium and high centering values for the damping/amplification parameter estimates for the 1-level of the distinguishing variable (e.g., values that indicate typical low, medium and high values of damping/amplification for the 1-level partner).
#' @param centFreqPart1 A vector of low, medium and high centering values for the partner influence on frequency parameter estimates for the 1-level of the distinguishing variable (e.g., values that indicate typical low, medium and high values of partner influence on frequency for the 1-level partner).
#' @param centDampPart1 A vector of low, medium and high centering values for the partner influence on damping/amplification parameter estimates for the 1-level of the distinguishing variable (e.g., values that indicate typical low, medium and high values of partner influence on damping/amplification for the 1-level partner).
#' @param sysVarName A name for the system variable being predicted (e.g., "Satisfaction").
#' @param dist A number indicating which level of the distinguishing variable to produce plots for. Must be either 0 or 1.
#' @param dist0name A name for the 0-level of the distinguishing variable (e.g., "Women").
#' @param dist1name A name for the 1-level of the distinguishing variable (e.g., "Men").

#' @import ggplot2
#' @export
coupleSysVarOutPlots <- function(basedata, centFreqSelf0, centDampSelf0, centFreqPart0, centDampPart0, centFreqSelf1, centDampSelf1, centFreqPart1, centDampPart1, sysVarName, dist, dist0name, dist1name)
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
	freqSData <- basedata
	freqSData$freqL <- basedata$obs_0 - centFreqSelf0[1]
	freqSData$freqH <- basedata$obs_0 - centFreqSelf0[3]
	freqSData$dampM <- basedata$d1_0 - centDampSelf0[2]
	freqSData$p_freqM <- basedata$p_obs_0 - centFreqPart0[2]
	freqSData$p_dampM <- basedata$p_d1_0 - centDampPart0[2]
	
	dampSData <- basedata
	dampSData$dampL <- basedata$d1_0 - centDampSelf0[1]
	dampSData$dampH <- basedata$d1_0 - centDampSelf0[3]
	dampSData$freqM <- basedata$obs_0 - centFreqSelf0[2]
	dampSData$p_freqM <- basedata$p_obs_0 - centFreqPart0[2]
	dampSData$p_dampM <- basedata$p_d1_0 - centDampPart0[2]
	
	freqPData <- basedata
	freqPData$p_freqL <- basedata$p_obs_0 - centFreqPart0[1]
	freqPData$p_freqH <- basedata$p_obs_0 - centFreqPart0[3]
	freqPData$freqM <- basedata$obs_0 - centFreqSelf0[2]
	freqPData$dampM <- basedata$d1_0 - centDampSelf0[2]
	freqPData$p_dampM <- basedata$p_d1_0 - centDampPart0[2]
	
	dampPData <- basedata
	dampPData$p_dampL <- basedata$p_d1_0 - centDampPart0[1]
	dampPData$p_dampH <- basedata$p_d1_0 - centDampPart0[3]
	dampPData$freqM <- basedata$obs_0 - centFreqSelf0[2]
	dampPData$dampM <- basedata$d1_0 - centDampSelf0[2]
	dampPData$p_freqM <- basedata$p_obs_0 - centFreqPart0[2]
	
	freqSL <- lm(sysVar ~ freqL + dampM + p_freqM + p_dampM, data=freqSData)
	freqSH <- lm(sysVar ~ freqH + dampM + p_freqM + p_dampM, data=freqSData)
	dampSL <- lm(sysVar ~ freqM + dampL + p_freqM + p_dampM, data=dampSData)
	dampSH <- lm(sysVar ~ freqM + dampH + p_freqM + p_dampM, data=dampSData)
	
	freqPL <- lm(sysVar ~ freqM + dampM + p_freqL + p_dampM, data=freqPData)
	freqPH <- lm(sysVar ~ freqM + dampM + p_freqH + p_dampM, data=freqPData)
	dampPL <- lm(sysVar ~ freqM + dampM + p_freqM + p_dampL, data=dampPData)
	dampPH <- lm(sysVar ~ freqM + dampM + p_freqM + p_dampH, data=dampPData)
		
	estFreqSL <- summary(freqSL)$coefficients[1,1]
	errFreqSL <- summary(freqSL)$coefficients[1,2]
	estFreqSH <- summary(freqSH)$coefficients[1,1]
	errFreqSH <- summary(freqSH)$coefficients[1,2]
	
	estDampSL <- summary(dampSL)$coefficients[1,1]
	errDampSL <- summary(dampSL)$coefficients[1,2]
	estDampSH <- summary(dampSH)$coefficients[1,1]
	errDampSH <- summary(dampSH)$coefficients[1,2]

	estFreqPL <- summary(freqPL)$coefficients[1,1]
	errFreqPL <- summary(freqPL)$coefficients[1,2]
	estFreqPH <- summary(freqPH)$coefficients[1,1]
	errFreqPH <- summary(freqPH)$coefficients[1,2]
	
	estDampPL <- summary(dampPL)$coefficients[1,1]
	errDampPL <- summary(dampPL)$coefficients[1,2]
	estDampPH <- summary(dampPH)$coefficients[1,1]
	errDampPH <- summary(dampPH)$coefficients[1,2]	
	
	freqSLdata <- c(1,1,estFreqSL,errFreqSL)
	freqSHdata <- c(1,2,estFreqSH,errFreqSH)
	dampSLdata <- c(2,1,estDampSL,errDampSL)
	dampSHdata <- c(2,2,estDampSH,errDampSH)
	
	freqPLdata <- c(3,1,estFreqPL,errFreqPL)
	freqPHdata <- c(3,2,estFreqPH,errFreqPH)
	dampPLdata <- c(4,1,estDampPL,errDampPL)
	dampPHdata <- c(4,2,estDampPH,errDampPH)
		
	plotData <- as.data.frame(rbind(freqSLdata, freqSHdata, dampSLdata, dampSHdata, freqPLdata, freqPHdata, dampPLdata, dampPHdata))
	names(plotData) <- c("Parameter","Level","sysVar","SE")
	plotData$Parameter <- factor(plotData$Parameter, labels=c("Freq-Self", "Damp-Self", "Freq-Part", "Damp-Part"))
	plotData$Level <- factor(plotData$Level, labels=c("Low","High"))
	plotData$errMin <- plotData[ ,3] - plotData[ ,4]
	plotData$errMax <- plotData[ ,3] + plotData[ ,4]
		
	Ymin <- min(basedata$sysVar, na.rm=T)
	Ymax <- max(basedata$sysVar, na.rm=T)
			
	couple <- ggplot(plotData, aes(x=Parameter, y=sysVar, fill=Level)) +
		geom_bar(position=position_dodge(), stat="identity") +
		geom_errorbar(aes(ymin=errMin, ymax=errMax), width=.2, position=position_dodge(.9)) +
		coord_cartesian(ylim=c(Ymin, Ymax)) +
		labs(title="Coupled Oscillator", subtitle= distName) +
		ylab(sysVarName) +
		theme(plot.title=element_text(size=12))
	
	} else if (dist==1)
	
	{
	freqSData <- basedata
	freqSData$freqL <- basedata$obs_1 - centFreqSelf1[1]
	freqSData$freqH <- basedata$obs_1 - centFreqSelf1[3]
	freqSData$dampM <- basedata$d1_1 - centDampSelf1[2]
	freqSData$p_freqM <- basedata$p_obs_1 - centFreqPart1[2]
	freqSData$p_dampM <- basedata$p_d1_1 - centDampPart1[2]
	
	dampSData <- basedata
	dampSData$dampL <- basedata$d1_1 - centDampSelf1[1]
	dampSData$dampH <- basedata$d1_1 - centDampSelf1[3]
	dampSData$freqM <- basedata$obs_1 - centFreqSelf1[2]
	dampSData$p_freqM <- basedata$p_obs_1 - centFreqPart1[2]
	dampSData$p_dampM <- basedata$p_d1_1 - centDampPart1[2]
	
	freqPData <- basedata
	freqPData$p_freqL <- basedata$p_obs_1 - centFreqPart1[1]
	freqPData$p_freqH <- basedata$p_obs_1 - centFreqPart1[3]
	freqPData$freqM <- basedata$obs_1 - centFreqSelf1[2]
	freqPData$dampM <- basedata$d1_1 - centDampSelf1[2]
	freqPData$p_dampM <- basedata$p_d1_1 - centDampPart1[2]
	
	dampPData <- basedata
	dampPData$p_dampL <- basedata$p_d1_1 - centDampPart1[1]
	dampPData$p_dampH <- basedata$p_d1_1 - centDampPart1[3]
	dampPData$freqM <- basedata$obs_1 - centFreqSelf1[2]
	dampPData$dampM <- basedata$d1_1 - centDampSelf1[2]
	dampPData$p_freqM <- basedata$p_obs_1 - centFreqPart1[2]
	
	freqSL <- lm(sysVar ~ freqL + dampM + p_freqM + p_dampM, data=freqSData)
	freqSH <- lm(sysVar ~ freqH + dampM + p_freqM + p_dampM, data=freqSData)
	dampSL <- lm(sysVar ~ freqM + dampL + p_freqM + p_dampM, data=dampSData)
	dampSH <- lm(sysVar ~ freqM + dampH + p_freqM + p_dampM, data=dampSData)
	
	freqPL <- lm(sysVar ~ freqM + dampM + p_freqL + p_dampM, data=freqPData)
	freqPH <- lm(sysVar ~ freqM + dampM + p_freqH + p_dampM, data=freqPData)
	dampPL <- lm(sysVar ~ freqM + dampM + p_freqM + p_dampL, data=dampPData)
	dampPH <- lm(sysVar ~ freqM + dampM + p_freqM + p_dampH, data=dampPData)
		
	estFreqSL <- summary(freqSL)$coefficients[1,1]
	errFreqSL <- summary(freqSL)$coefficients[1,2]
	estFreqSH <- summary(freqSH)$coefficients[1,1]
	errFreqSH <- summary(freqSH)$coefficients[1,2]
	
	estDampSL <- summary(dampSL)$coefficients[1,1]
	errDampSL <- summary(dampSL)$coefficients[1,2]
	estDampSH <- summary(dampSH)$coefficients[1,1]
	errDampSH <- summary(dampSH)$coefficients[1,2]

	estFreqPL <- summary(freqPL)$coefficients[1,1]
	errFreqPL <- summary(freqPL)$coefficients[1,2]
	estFreqPH <- summary(freqPH)$coefficients[1,1]
	errFreqPH <- summary(freqPH)$coefficients[1,2]
	
	estDampPL <- summary(dampPL)$coefficients[1,1]
	errDampPL <- summary(dampPL)$coefficients[1,2]
	estDampPH <- summary(dampPH)$coefficients[1,1]
	errDampPH <- summary(dampPH)$coefficients[1,2]	
	
	freqSLdata <- c(1,1,estFreqSL,errFreqSL)
	freqSHdata <- c(1,2,estFreqSH,errFreqSH)
	dampSLdata <- c(2,1,estDampSL,errDampSL)
	dampSHdata <- c(2,2,estDampSH,errDampSH)
	
	freqPLdata <- c(3,1,estFreqPL,errFreqPL)
	freqPHdata <- c(3,2,estFreqPH,errFreqPH)
	dampPLdata <- c(4,1,estDampPL,errDampPL)
	dampPHdata <- c(4,2,estDampPH,errDampPH)
		
	plotData <- as.data.frame(rbind(freqSLdata, freqSHdata, dampSLdata, dampSHdata, freqPLdata, freqPHdata, dampPLdata, dampPHdata))
	names(plotData) <- c("Parameter","Level","sysVar","SE")
	plotData$Parameter <- factor(plotData$Parameter, labels=c("Freq-Self", "Damp-Self", "Freq-Part", "Damp-Part"))
	plotData$Level <- factor(plotData$Level, labels=c("Low","High"))
	plotData$errMin <- plotData[ ,3] - plotData[ ,4]
	plotData$errMax <- plotData[ ,3] + plotData[ ,4]
		
	Ymin <- min(basedata$sysVar, na.rm=T)
	Ymax <- max(basedata$sysVar, na.rm=T)
			
	couple <- ggplot(plotData, aes(x=Parameter, y=sysVar, fill=Level)) +
		geom_bar(position=position_dodge(), stat="identity") +
		geom_errorbar(aes(ymin=errMin, ymax=errMax), width=.2, position=position_dodge(.9)) +
		coord_cartesian(ylim=c(Ymin, Ymax)) +
		labs(title="Coupled Oscillator", subtitle= distName) +
		ylab(sysVarName) +
		theme(plot.title=element_text(size=12))
	}
	plots <- couple
}


#' Plots the state variable's clo model predicted temporal trajectory 
#' 
#' Produces plots of the state variable's predicted temporal trajectory based on 3 versions of the clo model. The first is the "average" model, with all parameters centered on the sample averages. The other two are specified by the user by providing centering values for any desired parameters. For example, if prior analyses showed that the system variable was predicted by the damping parameter for the 0-level of the distinguishing variable, then a logical pair of models to compare would be one with that parameter centered at a low value, such as 1SD below the mean, and one with it centered at a high value, such as 1SD above the mean. All other parameters would be centered at the sample averages.
#' 
#' @param origdata A dataframe produced by the "dataPrep" function.
#' @param paramData A dataframe produced by the "indivClo" function.
#' @param paramM1 A vector of centering values for the first comparison model.
#' @param paramM2 A vector of centering values for the second comparison model.
#' @param dist0name A name for the 0-level of the distinguishing variable (e.g., "Women").
#' @param dist1name A name for the 1-level of the distinguishing variable (e.g., "Men").
#' @param obsName A name for the state variable (e.g., "Emotional Experience").
#' @param m1Name A name for the first user-specified model.
#' @param m2Name A name for the second user-specified model.

#' @import ggplot2
#' @export
cloPredTraj <- function(origdata, paramData, paramM1, paramM2, dist0name, dist1name, obsName, m1Name, m2Name)
{
	statedata0 <- origdata[origdata$dist0 == 1 & origdata$time ==1,] 
	start0 <- median(statedata0$obs_deTrend)
  	statedata1 <- origdata[origdata$dist0 == 0 & origdata$time ==1,] 
	start1 <- median(statedata1$obs_deTrend)
  	
  	maxtime <- quantile(origdata$time, prob=.75)
	plotTimes <- seq(1, maxtime, by=1)
	min <- min(origdata$obs_deTrend, na.rm=T)
	max <- max(origdata$obs_deTrend, na.rm=T)

	state <- c("y1"=start0, "y2"=0, "y3"=start1, "y4"=0)
	
	ave_obs_0 <- median(paramData$obs_0)
	ave_obs_1 <- median(paramData$obs_1)
	ave_p_obs_0 <- median(paramData$p_obs_0)
	ave_p_obs_1 <- median(paramData$p_obs_1)
	ave_d1_0 <- median(paramData$d1_0)
	ave_d1_1 <- median(paramData$d1_1)
	ave_p_d1_0 <- median(paramData$p_d1_0)
	ave_p_d1_1 <- median(paramData$p_d1_1)
	paramAve <- list(obs_0=ave_obs_0, obs_1=ave_obs_1, d1_0=ave_d1_0, d1_1=ave_d1_1, p_obs_0=ave_p_obs_0, p_obs_1=ave_p_obs_1, p_d1_0=ave_p_d1_0, p_d1_1=ave_p_d1_1)

	tempAve <- as.data.frame(deSolve::ode(y=state, times=plotTimes, func=cloCoupleOde, parms= paramAve))
	temp2 <- subset(tempAve, select=-c(y2, y4))
	names(temp2) <- c("time","d0.pred","d1.pred")
	temp3 <- reshape(temp2, direction='long', varying=c("d0.pred","d1.pred"), timevar="role", times=c("d0","d1"), v.names=c("pred"), idvar="time")
	temp3$roleNew <- factor(temp3$role, levels=c("d0","d1"), labels=c(dist0name, dist1name)) 
			
	plotData <- temp3[complete.cases(temp3), ]	
				
	aveCloPlot <- ggplot(plotData, aes(x=time)) +
				geom_line(aes(y= pred, color=roleNew), linetype="dotted", size= .8, na.rm=T) +
				scale_color_manual(name="Role", values=c("red","blue")) +
				ylab(obsName) +
				ylim(min, max) +
				labs(title="Predicted Trajectory", subtitle= "Average Model") +
				theme(plot.title=element_text(size=11))
	
	paramM1 <- ifelse(!is.na(paramM1), paramM1, paramAve)			
	tempM1 <- as.data.frame(deSolve::ode(y=state, times=plotTimes, func=cloCoupleOde, parms= paramM1))
	temp2 <- subset(tempM1, select=-c(y2, y4))
	names(temp2) <- c("time","d0.pred","d1.pred")
	temp3 <- reshape(temp2, direction='long', varying=c("d0.pred","d1.pred"), timevar="role", times=c("d0","d1"), v.names=c("pred"), idvar="time")
	temp3$roleNew <- factor(temp3$role, levels=c("d0","d1"), labels=c(dist0name, dist1name)) 
			
	plotData <- temp3[complete.cases(temp3), ]	
				
	m1CloPlot <- ggplot(plotData, aes(x=time)) +
				geom_line(aes(y= pred, color=roleNew), linetype="dotted", size= .8, na.rm=T) +
				scale_color_manual(name="Role", values=c("red","blue")) +
				ylab(obsName) +
				ylim(min, max) +
				labs(title="Predicted Trajectory", subtitle= m1Name) +
				theme(plot.title=element_text(size=11))

	paramM2 <- ifelse(!is.na(paramM2), paramM2, paramAve)	
	tempM2 <- as.data.frame(deSolve::ode(y=state, times=plotTimes, func=cloCoupleOde, parms= paramM2))
	temp2 <- subset(tempM2, select=-c(y2, y4))
	names(temp2) <- c("time","d0.pred","d1.pred")
	temp3 <- reshape(temp2, direction='long', varying=c("d0.pred","d1.pred"), timevar="role", times=c("d0","d1"), v.names=c("pred"), idvar="time")
	temp3$roleNew <- factor(temp3$role, levels=c("d0","d1"), labels=c(dist0name, dist1name)) 
			
	plotData <- temp3[complete.cases(temp3), ]	
				
	m2CloPlot <- ggplot(plotData, aes(x=time)) +
				geom_line(aes(y= pred, color=roleNew), linetype="dotted", size= .8, na.rm=T) +
				scale_color_manual(name="Role", values=c("red","blue")) +
				ylab(obsName) +
				ylim(min, max) +
				labs(title="Predicted Trajectory", subtitle= m2Name) +
				theme(plot.title=element_text(size=11))
					
	CLOplots <- list(aveCloPlot=aveCloPlot, m1CloPlot=m1CloPlot, m2CloPlot=m2CloPlot)
}




