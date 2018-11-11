
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
#' @param idConvention The number that was added to the dist0 partner to get the ID number for the dist1 partner.
#' @param dist0name A name for the level-0 of the distinguishing variable (e.g., "Women").
#' @param dist1name A name for the level-1 of the distinguishing variable (e.g., "Men").
#' @param obsName A name for the observed state variables being plotted (e.g., "Emotional Experience").
#' 
#' @return The function returns a list including: 1) the adjusted R^2 for the model for each dyad (called "R2"), 2) the parameter estimates for the model for each dyad (called "paramData", for use in either predicting, or being predicted by, the system variable), and 3) plots of the predicted values against the observed values for each dyad (called "plots"). The plots are also written to the working directory as a pdf file called "cloPlotsCouple.pdf"

#' @import ggplot2
#' @export
indivCloCouple <- function(basedata, idConvention, dist0name, dist1name, obsName)
{
	newDiD <- unique(factor(basedata$dyad))
	basedata <- basedata[complete.cases(basedata), ]
	min <- quantile(basedata$obs_deTrend, .1, na.rm=T)
	max <- quantile(basedata$obs_deTrend, .9,  na.rm=T)
	
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
			temp3$id <- ifelse(temp3$role == "d0", temp3$dyad, temp3$dyad + idConvention)
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

	results <- list(R2=r2, paramData=paramData, plots=plots)
}


#############################

#' Estimates an "un-coupled oscillator" model (e.g., only self frequency and damping with no coupling) for each dyad.
#' 
#' The second derivatives of the observed state variables (with linear trends removed) are predicted from each person's own observed state variable (again with linear trends removed), as well as each person's own first derivatives of the observed state variables (again with linear trends removed).
#'
#' @param basedata A dataframe that was produced with the "estDerivs" function.
#' @param idConvention The number that was added to the dist0 partner to get the ID number for the dist1 partner.
#' @param dist0name A name for the level-0 of the distinguishing variable (e.g., "Women").
#' @param dist1name A name for the level-1 of the distinguishing variable (e.g., "Men").
#' @param obsName A name for the observed state variables being plotted (e.g., "Emotional Experience").
#' 
#' @return The function returns a list including: 1) the adjusted R^2 for the model for each dyad (called "R2"), 2) the parameter estimates for the model for each dyad (called "paramData", for use in either predicting, or being predicted by, the system variable), and 3) plots of the predicted values against the observed values for each dyad (called "plots"). The plots are also written to the working directory as a pdf file called "cloPlotsUncouple.pdf"

#' @import ggplot2
#' @export
indivCloUncouple <- function(basedata, idConvention, dist0name, dist1name, obsName)
{
	newDiD <- unique(factor(basedata$dyad))
	basedata <- basedata[complete.cases(basedata), ]
	min <- quantile(basedata$obs_deTrend, .1, na.rm=T)
	max <- quantile(basedata$obs_deTrend, .9,  na.rm=T)
		
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
			temp3$id <- ifelse(temp3$role == "d0", temp3$dyad, temp3$dyad + idConvention)
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

	results <- list(R2=r2, paramData=paramData, plots=plots)
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
#' @return The function returns a named list including: 1) the adjusted R^2 for the uncoupled model for each dyad (called "R2uncouple"), 2) the adjusted R^2 for the coupled model for each dyad (called "R2couple"), and 3) the difference between the R-squares for each dyad (coupled - uncoupled, called "R2dif").

#' @export
cloIndivCompare <- function(basedata)
{
	newDiD <- unique(factor(basedata$dyad))

	R2uncouple <- vector()
	R2couple <- vector()
	R2dif <- vector()
	
		for (i in 1:length(newDiD))
		{
			datai <- basedata[basedata$dyad == newDiD[i], ]
			uncouple <- lm(d2 ~ dist0:obs_deTrend + dist0:d1 + dist1:obs_deTrend + dist1:d1 -1, na.action=na.exclude, data=datai)
			couple <- lm(d2 ~ dist0:obs_deTrend + dist0:d1 + dist0:p_obs_deTrend + dist0:p_d1 + dist1:obs_deTrend + dist1:d1 + dist1:p_obs_deTrend + dist1:p_d1-1, na.action=na.exclude, data=datai)	
			R2uncouple[[i]] <- summary(uncouple)$adj.r.squared
			R2couple[[i]] <- summary(couple)$adj.r.squared	
			R2dif[[i]] <- R2couple[[i]] - R2uncouple[[i]]
		}			
		output <- list(R2uncouple=R2uncouple, R2couple=R2couple, R2dif=R2dif)
		}


#' Compares a base line "intercept only" model to coupled and uncoupled oscillator models for predicting the system variable from the dynamic parameters.
#' 
#' The dynamic parameters used in these models come from the set of clo parameter: obs_0, obs_1, d1_0, d1_1, p_obs_0, p_obs_1, p_d1_0, p_d1_1. The 3 models compared are the baseline intercept-only, the uncoupled CLO (obs_0, obs_1, d1_0, d1_1), and the full coupled CLO. The system variable can be either dyadic (sysVarType = "dyad"), where both partners have the same score (e.g., relationship length) or individual (sysVarType = "indiv"), where the partners can have different scores (e.g., age). To make it easier to read the output, the dynamic parameters are renamed as follows: obs = freq, d1 = damp, p_obs = freqCoupling, p_d1 = dampCoupling.

#' 
#' @param basedata A dataframe containing the coupled oscillator parameter estimates produced by the "indivCloCouple" function.
#' @param sysVarType Whether the system variable is "dyad", which means both partners have the same socre, or "indiv" which means the partners can have different scores
#' @param dist0name A name for the level-0 of the distinguishing variable (e.g., "Women").
#' @param dist1name A name for the level-1 of the distinguishing variable (e.g., "Men").
#' @param sysVarName A name for the system variable being predicted (e.g., "Satisfaction").
#' 
#' @return The function returns a list including: 1) the lm or lme objects containing the full results for each model(called "models"), and 2) adjusted R^2 information for each model  (called "R2"). The function also displays histograms of the residuals and plots of the predicted values against observed values for each model. 

#' @export
cloSysVarOut <- function(basedata, sysVarType, dist0name, dist1name, sysVarName)
{
	basedata$dist1 <- ifelse(basedata$dist0 == 1, 0, 1)
	
	# Names for model parameters
	freq0name <- paste("freq",dist0name, sep="_")
	freq1name <- paste("freq",dist1name, sep="_")
	damp0name <- paste("damp",dist0name, sep="_")
	damp1name <- paste("damp",dist1name, sep="_")
	
	freqCouple0name <- paste("freqCoupling",dist0name, sep="_")
	freqCouple1name <- paste("freqCoupling",dist1name, sep="_")
	dampCouple0name <- paste("dampCoupling",dist0name, sep="_")
	dampCouple1name <- paste("dampCoupling",dist1name, sep="_")

	
	if(sysVarType != "indiv" & sysVarType != "dyad") 
    {
	stop("the sysVarType must be either indiv or dyad")
	}
	
	else if (sysVarType == "dyad")
 	{	
 	basedata <- basedata[!duplicated(basedata$dyad), ]

	base <- lm(sysVar ~ 1, data= basedata)
	names(base$coefficients) <- c("intercept")

	uncoupled <- lm(sysVar ~ obs_0 + d1_0 + obs_1 + d1_1, data= basedata)
	names(uncoupled$coefficients) <- c("intercept", freq0name, damp0name, freq1name, damp1name)
		
	coupled <- lm(sysVar ~ obs_0 + d1_0 + obs_1 + d1_1 + p_obs_0 + p_d1_0 + p_obs_1 + p_d1_1, data= basedata)
	names(coupled$coefficients) <- c("intercept", freq0name, damp0name, freq1name, damp1name, freqCouple0name, dampCouple0name, freqCouple1name, dampCouple1name)
		
	basePred <- predict(base)## all predictions will be the mean
	uncoupledPred <- predict(uncoupled)
	coupledPred <- predict(coupled)
	
	baseR2 <- summary(base)$adj.r.squared # will be zero
	uncoupledR2 <- summary(uncoupled)$adj.r.squared 
	coupledR2 <- summary(coupled)$adj.r.squared
	}
	
	else if (sysVarType == "indiv")
	{	
	base <- nlme::lme(sysVar ~ dist0, random= ~ 1 | dyad, data= basedata, na.action=na.omit, method="ML", control=nlme::lmeControl(opt="optim"))
	names(base$coefficients$fixed) <- c("intercept", dist0name)
	
	uncoupled <- nlme::lme(sysVar ~ dist0:obs_0 + dist0:d1_0 + dist1:obs_1 + dist1:d1_0 , random= ~ 1 | dyad, data= basedata, na.action=na.omit, method="ML", control=nlme::lmeControl(opt="optim"))
	names(uncoupled$coefficients$fixed) <- c("intercept", freq0name, damp0name, freq1name, damp1name)
		
	coupled <- nlme::lme(sysVar ~ dist0:obs_0 + dist0:d1_0 + dist1:obs_1 + dist1:d1_1 + 
								dist0:obs_1 + dist0:d1_1 + dist1:obs_0 + dist1:d1_0, random= ~ 1 | dyad, data= basedata, na.action=na.omit, method="ML", control=nlme::lmeControl(opt="optim"))
	names(coupled$coefficients$fixed) <- c("intercept", freq0name, damp0name, freq1name, damp1name, freqCouple0name, dampCouple0name, freqCouple1name, dampCouple1name)
	
	basePred <- predict(base)
	uncoupledPred <- predict(uncoupled)
	coupledPred <- predict(coupled)
	
	obs <- basedata$sysVar
	baseR2 <- summary(lm(obs ~ basePred))$adj.r.squared
	uncoupledR2 <- summary(lm(obs ~ uncoupledPred))$adj.r.squared
	coupledR2 <- summary(lm(obs ~ coupledPred))$adj.r.squared
	}
	
	
	min <- min(basedata$sysVar, na.rm=T)
	max <- max(basedata$sysVar, na.rm=T)
	
	ylabName <- paste(sysVarName, "observed", sep="_")
	xlabName <- paste(sysVarName, "predicted", sep="_")

	par(mfrow=c(2,1))
	hist(residuals(base))
	plot(basedata$sysVar ~ basePred, xlim=c(min, max), ylim=c(min, max), ylab=ylabName, xlab=xlabName, main="Baseline Model")

	par(mfrow=c(2,1))
	hist(residuals(uncoupled))
	plot(basedata$sysVar ~ uncoupledPred, xlim=c(min, max), ylim=c(min, max), ylab=ylabName, xlab=xlabName, main="Uncoupled Oscillator")
	
	par(mfrow=c(2,1))
	hist(residuals(coupled))
	plot(basedata$sysVar ~ coupledPred, xlim=c(min, max), ylim=c(min, max), ylab=ylabName, xlab=xlabName, main="Coupled Oscillator")
	

	models <- list(base=base, uncoupled=uncoupled, coupled=coupled)
	R2 <- list(baseR2=baseR2, uncoupledR2=uncoupledR2, coupledR2=coupledR2)
		
	output <- list(models=models, R2=R2)
}


#' Plots the state variable's clo model predicted temporal trajectory 
#' 
#' Produces plots of the state variable's predicted temporal trajectory based on 3 versions of the clo model. The first is the "average" model, with all parameters centered on the sample averages. The other two are specified by the user by providing centering values for any desired parameters. For example, if prior analyses showed that the system variable was predicted by the damping parameter for the level-0 of the distinguishing variable, then a logical pair of models to compare would be one with that parameter centered at a low value, such as 1SD below the mean, and one with it centered at a high value, such as 1SD above the mean. All other parameters would be centered at the sample averages.
#' 
#' @param origdata A dataframe produced by the "dataPrep" function.
#' @param paramData A dataframe produced by the "indivClo" function.
#' @param paramM1 A vector of centering values for the first comparison model.
#' @param paramM2 A vector of centering values for the second comparison model.
#' @param dist0name A name for the level-0 of the distinguishing variable (e.g., "Women").
#' @param dist1name A name for the level-1 of the distinguishing variable (e.g., "Men").
#' @param obsName A name for the state variable (e.g., "Emotional Experience").
#' @param m1Name A name for the first user-specified model.
#' @param m2Name A name for the second user-specified model.

#' @import ggplot2
#' @export

cloPredTraj <- function(origdata, paramData, paramM1, paramM2, dist0name, dist1name, obsName, m1Name, m2Name)
{
	statedata0 <- origdata[origdata$dist0 == 1 & origdata$time ==1,] 
	start0 <- median(statedata0$obs_deTrend, na.rm=T)
  	statedata1 <- origdata[origdata$dist0 == 0 & origdata$time ==1,] 
	start1 <- median(statedata1$obs_deTrend, na.rm=T)
  	
  	maxtime <- quantile(origdata$time, prob=.75)
	plotTimes <- seq(1, maxtime, by=1)
	min <- quantile(origdata$obs_deTrend, .1, na.rm=T)
	max <- quantile(origdata$obs_deTrend, .9,  na.rm=T)

	state <- c("y1"=start0, "y2"=0, "y3"=start1, "y4"=0)
	
	ave_obs_0 <- median(paramData$obs_0, na.rm=T)
	ave_obs_1 <- median(paramData$obs_1, na.rm=T)
	ave_p_obs_0 <- median(paramData$p_obs_0, na.rm=T)
	ave_p_obs_1 <- median(paramData$p_obs_1, na.rm=T)
	ave_d1_0 <- median(paramData$d1_0, na.rm=T)
	ave_d1_1 <- median(paramData$d1_1, na.rm=T)
	ave_p_d1_0 <- median(paramData$p_d1_0, na.rm=T)
	ave_p_d1_1 <- median(paramData$p_d1_1, na.rm=T)
	paramAve <- list(obs_0=ave_obs_0, obs_1=ave_obs_1, d1_0=ave_d1_0, d1_1=ave_d1_1, p_obs_0=ave_p_obs_0, p_obs_1=ave_p_obs_1, p_d1_0=ave_p_d1_0, p_d1_1=ave_p_d1_1)

	tempAve <- as.data.frame(deSolve::ode(y=state, times=plotTimes, func=cloCoupleOde, parms= paramAve))
	temp2 <- subset(tempAve, select=-c(y2, y4))
	names(temp2) <- c("time","d0.pred","d1.pred")
	temp3 <- reshape(temp2, direction='long', varying=c("d0.pred","d1.pred"), timevar="role", times=c("d0","d1"), v.names=c("pred"), idvar="time")
	temp3$roleNew <- factor(temp3$role, levels=c("d0","d1"), labels=c(dist0name, dist1name)) 
			
	plotData <- temp3[complete.cases(temp3), ]	
				
	aveCloPlot <- ggplot(plotData, aes(x=time)) +
				geom_line(aes(y= pred, color=roleNew), linetype="solid", size=1, na.rm=T) +
				scale_color_manual(name="Role", values=c("black","gray47")) +
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
				geom_line(aes(y= pred, color=roleNew), linetype="solid", size= 1, na.rm=T) +
				scale_color_manual(name="Role", values=c("black","gray47")) +
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
				geom_line(aes(y= pred, color=roleNew), linetype="solid", size= 1, na.rm=T) +
				scale_color_manual(name="Role", values=c("black","gray47")) +
				ylab(obsName) +
				ylim(min, max) +
				labs(title="Predicted Trajectory", subtitle= m2Name) +
				theme(plot.title=element_text(size=11))
					
	CLOplots <- list(aveCloPlot=aveCloPlot, m1CloPlot=m1CloPlot, m2CloPlot=m2CloPlot)
}


#' Compares a baseline "intercepts only" model to one including the system variable for predicting the dynamic parameters of the coupled oscillator model.
#' 
#' Multivariate correlated residuals dyadic models are used to predict the set of coupled oscillator parameters (obs_0, obs_1, d1_0, d1_1, p_obs_0, p_obs_1, p_d1_0, p_d1_1) from the system variable. The system variable can be either dyadic (sysVarType = "dyad"), where both partners have the same score (e.g., relationship length) or individual (sysVarType = "indiv"), where the partners can have different scores (e.g., age). If it is individual then both actor and partner effects of the system variable are included. To make it easier to read the output, the dynamic parameters are renamed as follows: obs = freq, d1 = damp, p_obs = freqCoupling, p_d1 = dampCoupling. For the individual model, a_ and p_ distinguish actor and partner effects respectively.
#' 
#' @param basedata A dataframe containing the coupled oscillator parameter estimates produced by the "indivCloCouple" function.
#' @param sysVarType Whether the system variable is "dyad", which means both partners have the same socre, or "indiv" which means the partners can have different scores
#' @param dist0name A name for the level-0 of the distinguishing variable (e.g., "Women").
#' @param dist1name A name for the level-1 of the distinguishing variable (e.g., "Men").
#' 
#' @return The function returns a list including: 1) the gls objects containing the full results for each model (called "models"), and 2) adjusted R^2 information for each model (called "R2").  

#' @export
cloSysVarIn <- function(basedata, sysVarType, dist0name, dist1name)
{
     basedata <- basedata[complete.cases(basedata), ] 
     
     # Names for intercepts
    freqInt0name <- paste("freq",dist0name, sep="_")
	freqInt1name <- paste("freq",dist1name, sep="_")
	dampInt0name <- paste("damp",dist0name, sep="_")
	dampInt1name <- paste("damp",dist1name, sep="_")
	freqCouplingInt0name <- paste("freqCoupling",dist0name, sep="_")
	freqCouplingInt1name <- paste("freqCoupling",dist1name, sep="_")
	dampCouplingInt0name <- paste("dampCoupling",dist0name, sep="_")
	dampCouplingInt1name <- paste("dampCoupling",dist1name, sep="_")

     
     # Names for indiv model parameters
	aFreq0name <- paste("a_freq_sysVar",dist0name, sep="_")
	pFreq0name <- paste("p_freq_sysVar",dist0name, sep="_")
	aFreq1name <- paste("a_freq_sysVar",dist1name, sep="_")
	pFreq1name <- paste("p_freq_sysVar",dist1name, sep="_")
	aDamp0name <- paste("a_damp_sysVar",dist0name, sep="_")
	pDamp0name <- paste("p_damp_sysVar",dist0name, sep="_")
	aDamp1name <- paste("a_damp_sysVar",dist1name, sep="_")
	pDamp1name <- paste("p_damp_sysVar",dist1name, sep="_")
	
	aFreqCoupling0name <- paste("a_freqCoupling_sysVar",dist0name, sep="_")
	pFreqCoupling0name <- paste("p_freqCoupling_sysVar",dist0name, sep="_")
	aFreqCoupling1name <- paste("a_freqCoupling_sysVar",dist1name, sep="_")
	pFreqCoupling1name <- paste("p_freqCoupling_sysVar",dist1name, sep="_")
	aDampCoupling0name <- paste("a_dampCoupling_sysVar",dist0name, sep="_")
	pDampCoupling0name <- paste("p_dampCoupling_sysVar",dist0name, sep="_")
	aDampCoupling1name <- paste("a_dampCoupling_sysVar",dist1name, sep="_")
	pDampCoupling1name <- paste("p_dampCoupling_sysVar",dist1name, sep="_")


	# Names for dyad model parameters

	freqSysVar0name <- paste("freq_SysVar",dist0name, sep="_")
	freqSysVar1name <- paste("freq_SysVar",dist1name, sep="_")
	dampSysVar0name <- paste("damp_SysVar",dist0name, sep="_")
	dampSysVar1name <- paste("damp_SysVar",dist1name, sep="_")
	freqCouplingSysVar0name <- paste("freqCoupling_SysVar",dist0name, sep="_")
	freqCouplingSysVar1name <- paste("freqCoupling_SysVar",dist1name, sep="_")
	dampCouplingSysVar0name <- paste("dampCoupling_SysVar",dist0name, sep="_")
	dampCouplingSysVar1name <- paste("dampCoupling_SysVar",dist1name, sep="_")

      if(sysVarType != "indiv" & sysVarType != "dyad") 
   {
			stop("the sysVarType must be either indiv or dyad")
	}
 	else if(sysVarType=="indiv")
   {
	## Case with individual level sysVar

	data1 <- subset(basedata, select=c(dyad, sysVar, dist0))
	data2 <- stats::reshape(data1, idvar="dyad", timevar = "dist0", direction= "wide")
	data3 <- subset(basedata, select=c(dyad, obs_0:p_d1_1))
	data4 <- data3[!duplicated(data3$dyad), ]
	data5 <- plyr::join(data4, data2)
	colnames(data5) <- c("dyad", "paramEst1", "paramEst2", "paramEst3", "paramEst4", "paramEst5", "paramEst6", "paramEst7", "paramEst8","sysVar1", "sysVar0")

	data6 <- stats::reshape(data5, varying=c("paramEst1", "paramEst2", "paramEst3", "paramEst4","paramEst5", "paramEst6", "paramEst7", "paramEst8"), timevar="parameter", idvar="dyad", direction="long", sep="")
	data7 <- data6[complete.cases(data6), ] 

	data7$parameter <- factor(data7$parameter, levels = c(1:8), labels=c("freq0","damp0","freqCoupling0", "dampCoupling0", "freq1","damp1","freqCoupling1", "dampCoupling1"))
	data7$freq0 <- ifelse(data7$parameter == "freq0", 1, 0) 
	data7$damp0 <- ifelse(data7$parameter == "damp0", 1, 0) 
	data7$freqCoupling0 <- ifelse(data7$parameter == "freqCoupling0", 1, 0) 
	data7$dampCoupling0 <- ifelse(data7$parameter == "dampCoupling0", 1, 0) 
	data7$freq1 <- ifelse(data7$parameter == "freq1", 1, 0) 
	data7$damp1 <- ifelse(data7$parameter == "damp1", 1, 0) 
	data7$freqCoupling1 <- ifelse(data7$parameter == "freqCoupling1", 1, 0) 
	data7$dampCoupling1 <- ifelse(data7$parameter == "dampCoupling1", 1, 0) 


		base <- nlme::gls(paramEst ~ freq0 + damp0 + freqCoupling0 + dampCoupling0 + freq1 + damp1 + freqCoupling1 + dampCoupling1 -1, correlation=nlme::corCompSymm(form= ~ 1 |dyad), weights=nlme::varIdent(form = ~ 1|parameter), na.action=na.omit, method="ML", control=nlme::lmeControl(opt="optim"), data=data7)
	names(base$coefficients) <- c(freqInt0name, dampInt0name, freqCouplingInt0name, dampCouplingInt0name, freqInt1name, dampInt1name, freqCouplingInt1name, dampCouplingInt1name)


	sysVarIn <- nlme::gls(paramEst ~ freq0 + damp0 + freqCoupling0 + dampCoupling0 + freq1 + damp1 + freqCoupling1 + dampCoupling1 +
		freq0:sysVar0 + freq0:sysVar1 + damp0:sysVar0 + damp0:sysVar1 + freqCoupling0:sysVar0 + freqCoupling0:sysVar1 + dampCoupling0:sysVar0 + dampCoupling0:sysVar1 + 
	freq1:sysVar1 + freq1:sysVar0 + damp1:sysVar1 + damp1:sysVar0 + freqCoupling1:sysVar1 + freqCoupling1:sysVar0 + dampCoupling1:sysVar1 + dampCoupling1:sysVar0 -1,
	correlation=nlme::corCompSymm(form= ~ 1 |dyad), weights=nlme::varIdent(form = ~ 1|parameter), na.action=na.omit, method="ML", control=nlme::lmeControl(opt="optim"), data=data7)
	names(sysVarIn$coefficients) <- c(freqInt0name, dampInt0name, freqCouplingInt0name, dampCouplingInt0name, freqInt1name, dampInt1name, freqCouplingInt1name, dampCouplingInt1name,   aFreq0name, pFreq0name, aDamp0name, pDamp0name, aFreqCoupling0name, pFreqCoupling0name, aDampCoupling0name, pDampCoupling0name,   aFreq1name, pFreq1name, aDamp1name, pDamp1name, aFreqCoupling1name, pFreqCoupling1name, aDampCoupling1name, pDampCoupling1name )

	
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
	data2 <- subset(data1, select=-c(id, dist0))
	colnames(data2) <- c("paramEst1", "paramEst2", "paramEst3", "paramEst4","paramEst5","paramEst6","paramEst7","paramEst8", "dyad", "sysVar")
	data3 <- stats::reshape(data2, varying=c("paramEst1", "paramEst2", "paramEst3", "paramEst4","paramEst5","paramEst6","paramEst7","paramEst8"), timevar="parameter", idvar="dyad", direction="long", sep="")

	data3$parameter <- factor(data3$parameter, levels = c(1:8), labels=c("freq0","damp0","freqCoupling0", "dampCoupling0", "freq1","damp1","freqCoupling1", "dampCoupling1"))
	data3$freq0 <- ifelse(data3$parameter == "freq0", 1, 0) 
	data3$damp0 <- ifelse(data3$parameter == "damp0", 1, 0) 
	data3$freqCoupling0 <- ifelse(data3$parameter == "freqCoupling0", 1, 0) 
	data3$dampCoupling0 <- ifelse(data3$parameter == "dampCoupling0", 1, 0) 
	data3$freq1 <- ifelse(data3$parameter == "freq1", 1, 0) 
	data3$damp1 <- ifelse(data3$parameter == "damp1", 1, 0) 
	data3$freqCoupling1 <- ifelse(data3$parameter == "freqCoupling1", 1, 0) 
	data3$dampCoupling1 <- ifelse(data3$parameter == "dampCoupling1", 1, 0) 


	base <- nlme::gls(paramEst ~ freq0 + damp0 + freqCoupling0 + dampCoupling0 + freq1 + damp1 + freqCoupling1 + dampCoupling1 -1, correlation=nlme::corCompSymm(form= ~ 1 |dyad), weights=nlme::varIdent(form = ~ 1|parameter), na.action=na.omit, method="ML", control=nlme::lmeControl(opt="optim"), data=data3)
	names(base$coefficients) <- c(freqInt0name, dampInt0name, freqCouplingInt0name, dampCouplingInt0name, freqInt1name, dampInt1name, freqCouplingInt1name, dampCouplingInt1name)

	sysVarIn <- nlme::gls(paramEst ~ freq0 + damp0 + freqCoupling0 + dampCoupling0 + freq1 + damp1 + freqCoupling1 + dampCoupling1 + sysVar:freq0 + sysVar:damp0 + sysVar:freqCoupling0 + sysVar:dampCoupling0 + sysVar:freq1 + sysVar:damp1 + sysVar:freqCoupling1 + sysVar:dampCoupling1 -1, correlation=nlme::corCompSymm(form= ~ 1 |dyad), weights=nlme::varIdent(form = ~ 1|parameter), na.action=na.omit, method="ML", control=nlme::lmeControl(opt="optim"), data=data3)
	names(sysVarIn$coefficients) <- c(freqInt0name, dampInt0name, freqCouplingInt0name, dampCouplingInt0name, freqInt1name, dampInt1name, freqCouplingInt1name, dampCouplingInt1name, freqSysVar0name, dampSysVar0name, freqCouplingSysVar0name, dampCouplingSysVar0name, freqSysVar1name, dampSysVar1name, freqCouplingSysVar1name, dampCouplingSysVar1name)
		
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

