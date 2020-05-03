
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

############## estDerivs

#' Estimates first and second derivatives of an oberved state variable
#'
#' This function makes use of 2 functions written by Steven Boker, "gllaWMatrix" and "gllaEmbed" which are available on his website, http://people.virginia.edu/~smb3u/. It fits a coupled oscillator model for each dyad at different combinations of the input parameters (tau, embeds) and returns the input values and period of oscillation that maximize the R^2 for each dyad. It also estimates first and second derivatives of the observed state variable for each person at the input values that maximize the R^2 for that dyad and returns a dataframe that contains them.
#'
#' @param prepData A dataframe that was produced with the "dataPrep" function.
#' @param taus A vector containing the values of tau to use. Tau indicates the number of time points to lag in the lagged data matrix (see Boker, S.M., Deboeck, P.R., Edler, C., & Keel, P.K. (2010). Generalized local linear approximation of derivatives from time series. In S.M. Chow & E. Ferrer (Eds.), Statistical Methods for Modeling Human Dynamics: An Interdisciplinary Dialogue (pp. 161-178). New York, NY: Taylor & Francis Group). The first derivative is estimated as the mean of the two adjacent slopes across that number of lags, e.g., if tau = 2 then the estimate of the first derivative at time = t is based on the mean of the slopes left and right of time t across 2 observations each. The second derivative is the difference in the two slopes with respect to time. Tau = 1 is sensitive to noise and increasing its value acts as smoothing. 
#' @param embeds A vector containing the values of embeds to use. Embeds indicates the number of columns in the lagged data matrix. The minimum = 3 for 2nd order derivatives and higher values increase smoothing.
#' @param delta A value indicating the inter-observation interval. For example, if delta = 2, then every second observation is used in the estimation process.
#' @param idConvention The value that was added to the dist1 ID number to get the dist2 ID number
#' @examples
#' data <- rties_ExampleDataShort
#' newData <- dataPrep(basedata=data, dyadId="couple", personId="person", 
#' obs_name="dial", dist_name="female", time_name="time")
#' taus <-c(2,3)
#' embeds <- c(3,4)
#' delta <- 1
#' derivs <- estDerivs(prepData=newData, taus=taus, embeds=embeds, delta=delta, idConvention=500)
#' head(derivs$fitTable)
#' summary(derivs$fitTable[ ,4]) # summary of R-square
#' summary(derivs$fitTable[ ,5]) # summary of period of oscillation
#' 
#' @return The function returns a list including: 1) "data" which is a dataframe containing first and second derivative estimates of an observed state variable, and 2) "fitTable" which shows the maximal R^2 achieved for each dyad for a coupled oscillator model, along with the associated tau, embed and estimated period of oscillation.

#' @export

estDerivs <- function(prepData, taus, embeds, delta, idConvention)
{   
  basedata <- prepData
  basedata <- basedata[stats::complete.cases(basedata), ] 
  params <- expand.grid(taus=taus, embeds=embeds)
  dyadId <- unique(factor(basedata$dyad))

  derivData <- list()
  fitTable <- list()

  for(d in 1:length(dyadId)){
    r <- list()
    freq0 <- list()
    freq1 <- list()
	dataiFull <- basedata[basedata$dyad == dyadId[d],]
	datai <- dataiFull[ which(dataiFull$dist0 == 1), ]

    for(i in 1:nrow(params)){
	  # Estimate derivatives for each parameter combination
	  obsEmbed <- gllaEmbed(datai$obs_deTrend, tau=params[i,1], embed=params[i,2])
	  p_obsEmbed <- gllaEmbed(datai$p_obs_deTrend, tau=params[i,1], embed=params[i,2])   

	  obsLLA <- gllaWMatrix(tau=params[i,1], embed=params[i,2], deltaT=delta, order=2)
	  obsDeriv <- obsEmbed[,1:dim(obsEmbed)[2]] %*% obsLLA
	  p_obsDeriv <- p_obsEmbed[,1:dim(p_obsEmbed)[2]] %*% obsLLA

	  idLength <- dim(obsDeriv)[1]
	  dist0 <- rep(unique(datai$dist0), idLength)
	  dist1 <- rep(unique(datai$dist1), idLength)

	  deriv0 <- cbind(obsDeriv, dist0)
	  deriv1 <- cbind(p_obsDeriv, dist1)
	  deriv0full <- cbind(deriv0, deriv1)
	  dimnames(deriv0full) <- list(NULL, c("obs_deTrend","d1","d2","dist0","p_obs_deTrend","p_d1","p_d2","dist1")) 
	  deriv1full <- cbind(deriv1, deriv0)
	  dimnames(deriv1full) <- list(NULL, c("obs_deTrend","d1","d2","dist0","p_obs_deTrend","p_d1","p_d2","dist1"))

	  derivi <- rbind(deriv1full, deriv0full)
	  derivi <- as.data.frame(derivi)

	  # fit CLO and get R^2 for that combination of parameters
	  out <- stats::lm(d2 ~ dist0:obs_deTrend + dist0:d1 + dist0:p_obs_deTrend + dist0:p_d1 + dist1:obs_deTrend + dist1:d1 + dist1:p_obs_deTrend + dist1:p_d1 -1, data=derivi)
	  r[i] <- summary(out)$adj.r.squared 
	  freq0[i] <- out$coefficients[1]
	  freq1[i] <- out$coefficients[5]
	 }
	  
	## get maximum R^2 for given dyad and the corresponding tau & embed values
	maxR <- max(unlist(r))
	paramRow <- which(r==maxR)
	select <- params[paramRow,]
	tau <- select[1,1]
	embed <- select[1,2]
	temp0 <- unlist(freq0)
	freq0 <- temp0[paramRow]
	n0 <- abs(as.numeric(freq0))
	temp1 <- unlist(freq1)
	freq1 <- temp1[paramRow]
	n1 <- abs(as.numeric(freq1))

	## get period (1 cycle (peak to peak)) in given time units associated with highest R^2

	if (freq0 >= 0 | freq1 >= 0) {print("error: frequency parameter is not negative")}
	  
	period0 <- (2*pi) / (sqrt(n0))
	period1 <- (2*pi) / (sqrt(n1))

    # Estimate derivatives with selected tau and embed for each dyad
	obsEmbed <- gllaEmbed(datai$obs_deTrend, tau=tau, embed=embed)
	p_obsEmbed <- gllaEmbed(datai$p_obs_deTrend, tau=tau, embed=embed)   

	obsLLA <- gllaWMatrix(tau=tau, embed=embed, deltaT=delta, order=2)
	obsDeriv <- obsEmbed[,1:dim(obsEmbed)[2]] %*% obsLLA
	p_obsDeriv <- p_obsEmbed[,1:dim(p_obsEmbed)[2]] %*% obsLLA
	   
	idLength <- dim(obsDeriv)[1]
	dist0 <- rep(unique(datai$dist0), idLength)
	dist1 <- rep(unique(datai$dist1), idLength)
	time <- seq_len(idLength)
	dyad <- rep(unique(datai$dyad, idLength))
	deriv0 <- cbind(dyad, dist0, time, obsDeriv)
	deriv1 <- cbind(dyad, dist1, time, p_obsDeriv)

    deriv0full <- cbind(deriv0, deriv1)
    id0 <- rep(unique(datai$dyad + idConvention), idLength)
    deriv0full <- cbind(deriv0full, id0)
    dimnames(deriv0full) <- list(NULL, c("dyad","dist0", "time","obs_deTrend","d1","d2","p_dyad","dist1","p_time","p_obs_deTrend","p_d1","p_d2", "id"))	  
	
	deriv1full <- cbind(deriv1, deriv0)	
	id1 <- rep(unique(datai$dyad), idLength)
	deriv1full <- cbind(deriv1full, id1)
    dimnames(deriv1full) <- list(NULL, c("dyad","dist0", "time","obs_deTrend","d1","d2","p_dyad","dist1","p_time","p_obs_deTrend","p_d1","p_d2", "id"))	  

    deriv <- rbind(deriv1full, deriv0full)
	deriv <- as.data.frame(deriv)   
	derivData[[d]] <- deriv
	fitTable[[d]] <- c("dyad"= unique(deriv$dyad), "tau"=tau, "embed"= embed, 
	   "Rsqr"= maxR, "Period0"=period0, "Period1"=period1)
  }

   ## output, which includes the derivative data and the fit table

  derivD <- as.data.frame(do.call(rbind, derivData))
  fitTable <- as.data.frame(do.call(rbind, fitTable))
  fitTable <- round(fitTable, 2)

  derivOut <- list("data"=derivD, "fitTable"= fitTable)
}

################ cloCoupleOde

#' Provides the equation for a coupled oscillator model for the differential equation solver (ode) to plot

#' @param t A parameter used by the ode and passed by functions calling cloCoupleOde
#' @param state Another parameter used by the ode and passed by functions calling cloCoupleOde
#' @param parameters Another parameter used by the ode and passed by functions calling cloCoupleOde
#' 
#' @return A list with the rates of change for each state variable.

cloCoupledOde <- function(t, state, parameters)
{
  with(as.list(c(state, parameters)), {
	dy1 <- y2
	dy2 <- y1*obs_0 + y2*d1_0 + y3*p_obs_0 + y4*p_d1_0
	dy3 <- y4
	dy4 <- y3*obs_1 + y4*d1_1 + y1*p_obs_1 + y2*p_d1_1
	list(c(dy1, dy2, dy3, dy4))		
  })
}

############### cloUncoupledOde

#' Provides the equation for an un-coupled oscillator model for the differential equation solver (ode) to plot

#' @param t A parameter used by the ode and passed by functions calling cloCoupleOde
#' @param state Another parameter used by the ode and passed by functions calling cloCoupleOde
#' @param parameters Another parameter used by the ode and passed by functions calling cloCoupleOde
#' #' 
#' @return A list with the rates of change for each state variable.

cloUncoupledOde <- function(t, state, parameters)
{
  with(as.list(c(state, parameters)), {
	dy1 <- y2
	dy2 <- y1*obs_0 + y2*d1_0 
	dy3 <- y4
	dy4 <- y3*obs_1 + y4*d1_1 
	list(c(dy1, dy2, dy3, dy4))		
  })
}

################### indivClo

#' Estimates either an uncoupled or coupled oscillator model for each dyad.
#' 
#' Both models predict the second derivatives of the observed state variables (with linear trends removed). For the uncoupled oscillator, the predictors are each person's own observed state variables (again with linear trends removed), as well as each person's own first derivatives of the observed state variables (again with linear trends removed. For the coupled oscillator, the predictors are each person's own and partner's observed state variables (again with linear trends removed), as well as each person's own and partner's first derivatives of the observed state variables (again with linear trends removed).
#'
#' @param derivData A dataframe that was produced with the "estDerivs" function.
#' @param whichModel Whether the model to be estimated is the "uncoupled" or "coupled" oscillator.
#' @examples
#' data <- rties_ExampleDataShort
#' newData <- dataPrep(basedata=data, dyadId="couple", personId="person", 
#' obs_name="dial", dist_name="female", time_name="time")
#' taus <-c(2,3)
#' embeds <- c(3,4)
#' delta <- 1
#' derivs <- estDerivs(prepData=newData, taus=taus, embeds=embeds, delta=delta, idConvention=500)
#' clo <- indivClo(derivData=derivs$data, whichModel="coupled")
#' summary(clo$R2)
#' head(clo$params)

#' @return The function returns a list including: 1) the adjusted R^2 for the model for each dyad (called "R2"), and 2) the parameter estimates for the model for each dyad (called "params", for use in either predicting, or being predicted by, the system variable).

#' @export

indivClo <- function(derivData, whichModel)
{
  basedata <- derivData
  
  param <- list()
  
  if(whichModel != "uncoupled" & whichModel != "coupled") {
  	stop("the model type must be either uncoupled or coupled")
	
	} else if (whichModel == "uncoupled"){
	  model <- stats::formula(d2 ~ dist0:obs_deTrend + dist0:d1 + dist1:obs_deTrend + dist1:d1 -1)
	  obs_0 <- param[1]	
	  d1_0 <- param[2]
	  obs_1 <- param[3]
	  d1_1 <- param[4]
	  paramClo <- c("obs_0"= obs_0, "d1_0"= d1_0, "obs_1"=obs_1, "d1_1"=d1_1)
	  paramNames <- c("obs_0","d1_0","obs_1","d1_1","dyad")

      } else if (whichModel == "coupled"){
      	model <- stats::formula(d2 ~ dist0:obs_deTrend + dist0:d1 + dist0:p_obs_deTrend + dist0:p_d1 + dist1:obs_deTrend + dist1:d1 + dist1:p_obs_deTrend + dist1:p_d1 -1)
      	obs_0 <- param[1]	
		d1_0 <- param[2]
		p_obs_0 <- param[3]
		p_d1_0 <- param[4]
	    obs_1 <- param[5]
		d1_1 <- param[6]
		p_obs_1 <- param[7]
		p_d1_1 <- param[8]
		paramClo <- c("obs_0"= obs_0, "d1_0"= d1_0, "p_obs_0"= p_obs_0, "p_d1_0"=p_d1_0, "obs_1"=obs_1, "d1_1"=d1_1, "p_obs_1"= p_obs_1, "p_d1_1"= p_d1_1)
		paramNames <- c("obs_0","d1_0","p_obs_0","p_d1_0","obs_1","d1_1","p_obs_1","p_d1_1","dyad")
  }	

  newDiD <- unique(factor(basedata$dyad))
  basedata <- basedata[stats::complete.cases(basedata), ]
  R2 <- vector()
  	
  for (i in 1:length(newDiD)){
    datai <- basedata[basedata$dyad == newDiD[i], ]
	m <- stats::lm(model, na.action=na.exclude, data=datai)
	R2[[i]] <- summary(m)$adj.r.squared
	param[[i]] <- round(as.numeric(m$coefficients), 5)
	numParam <- length(m$coefficients)
	param[[i]][numParam + 1] <- unique(datai$dyad)
  }			
  params <- as.data.frame(do.call(rbind, param))
  colnames(params) <- paramNames	
  results <- list(R2=R2, params=params)
}


################### indivCloCompare

#' Compares model fit for the uncoupled and coupled oscillator for each dyad's state trajectories using an R-square comparison. 
#' 
#' Fits an uncoupled and coupled oscillator model to each dyad's observed state variables and returns the adjusted R-squares, along with the difference between them (coupled - uncoupled, so positive values indicate better fit for the more complex model).
#'
#' @param derivData A dataframe that was produced with the "estDerivs" function.
#' #' @examples
#' data <- rties_ExampleDataShort
#' newData <- dataPrep(basedata=data, dyadId="couple", personId="person", 
#' obs_name="dial", dist_name="female", time_name="time")
#' taus <-c(2,3)
#' embeds <- c(3,4)
#' delta <- 1
#' derivs <- estDerivs(prepData=newData, taus=taus, embeds=embeds, delta=delta, idConvention=500)
#' compare <- indivCloCompare(derivData=derivs$data)
#' summary(compare$R2couple)
#' 
#' @return The function returns a named list including: 1) the adjusted R^2 for the uncoupled model for each dyad (called "R2uncouple"), 2) the adjusted R^2 for the coupled model for each dyad (called "R2couple"), and 3) the difference between the R-squares for each dyad (coupled - uncoupled, called "R2dif").

#' @export

indivCloCompare <- function(derivData)
{
  basedata <- derivData
  
  newDiD <- unique(factor(basedata$dyad))
  R2uncouple <- vector()
  R2couple <- vector()
  R2dif <- vector()
  
  for (i in 1:length(newDiD)){
	datai <- basedata[basedata$dyad == newDiD[i], ]
	m1 <- stats::formula(d2 ~ dist0:obs_deTrend + dist0:d1 + dist1:obs_deTrend + dist1:d1 -1)
	uncouple <- stats::lm(m1, na.action=na.exclude, data=datai)
	m2 <- stats::formula(d2 ~ dist0:obs_deTrend + dist0:d1 + dist0:p_obs_deTrend + dist0:p_d1 + dist1:obs_deTrend + dist1:d1 + dist1:p_obs_deTrend + dist1:p_d1-1)
	couple <- stats::lm(m2, na.action=na.exclude, data=datai)	
	R2uncouple[[i]] <- summary(uncouple)$adj.r.squared
	R2couple[[i]] <- summary(couple)$adj.r.squared	
	R2dif[[i]] <- R2couple[[i]] - R2uncouple[[i]]
  }			
		
  output <- list(R2uncouple=R2uncouple, R2couple=R2couple, R2dif=R2dif)
}


################ indivCloPlots

#' Produces plots of either an uncoupled or coupled oscillator model-predicted trajectories overlaid on raw data for each dyad.
#' 
#' The observed and CLO-model predicted state variables (with linear trends removed) are plotted for each dyad individually.  
#'
#' @param derivData A dataframe that was produced with the "estDerivs" function.
#' @param whichModel Whether the model to be estimated is the "uncoupled" or "coupled" oscillator.
#' @param idConvention The number that was added to the dist0 partner to get the ID number for the dist1 partner.
#' @param dist0name An optional name for the level-0 of the distinguishing variable (e.g., "Women"). Default is dist0.
#' @param dist1name An optional name for the level-1 of the distinguishing variable (e.g., "Men"). Default is dist1.
#' @param plot_obs_name An optional name for the observed state variables being plotted (e.g., "Emotional Experience"). Default is observed.
#' @param minMax An optional vector with desired minimum and maximum quantiles to be used for setting the y-axis range on the plots, e.g., minMax <- c(.1, .9) would set the y-axis limits to the 10th and 90th percentiles of the observed state variables. Default is to use the minimum and maximum observed values of the state variables.
#' @param printPlots If true (the default) plots are displayed on the screen.
#' @examples
#' # See vignettes for examples.
#' 
#' @return A list plots of the predicted values against the observed values for each dyad.

#' @import ggplot2
#' @export

indivCloPlots <- function(derivData, whichModel, idConvention, dist0name=NULL, dist1name=NULL, plot_obs_name=NULL, minMax=NULL, printPlots=T)
{
  basedata <- derivData
  
  param <- list()
  
  if(is.null(dist0name)){dist0name <- "dist0"}
  if(is.null(dist1name)){dist1name <- "dist1"}
  if(is.null(plot_obs_name)){plot_obs_name <- "observed"}
  
  if(is.null(minMax)){
    min <- min(basedata$obs_deTrend, na.rm=T)
    max <- max(basedata$obs_deTrend, na.rm=T)
  } else {
    min <- stats::quantile(basedata$obs_deTrend, minMax[1], na.rm=T)
    max <- stats::quantile(basedata$obs_deTrend, minMax[2],  na.rm=T)
  }
  
  if(whichModel != "uncoupled" & whichModel != "coupled") {
  	stop("the model type must be either uncoupled or coupled")
	
	} else if (whichModel == "uncoupled"){
	  model <- stats::formula(d2 ~ dist0:obs_deTrend + dist0:d1 + dist1:obs_deTrend + dist1:d1 -1)
	  obs_0 <- param[1]	
	  d1_0 <- param[2]
	  obs_1 <- param[3]
	  d1_1 <- param[4]
	  paramClo <- c("obs_0"= obs_0, "d1_0"= d1_0, "obs_1"=obs_1, "d1_1"=d1_1)
	  paramNames <- c("obs_0","d1_0","obs_1","d1_1","dyad")
	  odeFunction <- cloUncoupledOde

      } else if (whichModel == "coupled"){
      	model <- stats::formula(d2 ~ dist0:obs_deTrend + dist0:d1 + dist0:p_obs_deTrend + dist0:p_d1 + dist1:obs_deTrend + dist1:d1 + dist1:p_obs_deTrend + dist1:p_d1 -1)
      	obs_0 <- param[1]	
		d1_0 <- param[2]
		p_obs_0 <- param[3]
		p_d1_0 <- param[4]
	    obs_1 <- param[5]
		d1_1 <- param[6]
		p_obs_1 <- param[7]
		p_d1_1 <- param[8]
		paramClo <- c("obs_0"= obs_0, "d1_0"= d1_0, "p_obs_0"= p_obs_0, "p_d1_0"=p_d1_0, "obs_1"=obs_1, "d1_1"=d1_1, "p_obs_1"= p_obs_1, "p_d1_1"= p_d1_1)
		paramNames <- c("obs_0","d1_0","p_obs_0","p_d1_0","obs_1","d1_1","p_obs_1","p_d1_1","dyad")
		odeFunction <- cloCoupledOde
  }	
   
  newDiD <- unique(factor(basedata$dyad))
  basedata <- basedata[stats::complete.cases(basedata), ]
  plots <- list()
	
  for (i in 1:length(newDiD)){
    statedatai <- basedata[basedata$dyad == newDiD[i] & basedata$dist0 == 0,] 
  	maxtime <- max(statedatai$time) 
 	plotTimes <- seq(1, maxtime, by=1)
 	time <- obs_deTrend <- p_obs_deTrend <- NULL
 	start <- suppressWarnings(subset(statedatai, time==c(1:5), select=c(obs_deTrend, p_obs_deTrend)))
	y1 <- mean(start$obs_deTrend, na.rm=T)
 	y2 <- 0
	y3 <- mean(start$p_obs_deTrend, na.rm=T)
 	y4 <- 0
 	statei <- c("y1"=y1, "y2"=y2, "y3"=y3, "y4"=y4)
			
	datai <- basedata[basedata$dyad == newDiD[i], ]
	m <- stats::lm(model, na.action=na.exclude, data=datai)
	param[[i]] <- round(as.numeric(m$coefficients), 5)
	numParam <- length(m$coefficients)
	param[[i]][numParam + 1] <- unique(datai$dyad)
    names(param[[i]]) <- paramNames

	temp <- as.data.frame(deSolve::ode(y=statei, times=plotTimes, func= odeFunction, parms= param[[i]]))
	vars1 <- c("y2", "y4")
	temp2 <- temp[ ,!(names(temp) %in% vars1)]
	names(temp2) <- c("time","d0.pred","d1.pred")
	temp2$dyad <- statedatai$dyad
	temp3 <- stats::reshape(temp2, direction='long', varying=c("d0.pred","d1.pred"), timevar="role", times=c("d0","d1"), v.names=c("pred"), idvar="time")
	temp3$id <- ifelse(temp3$role == "d0", temp3$dyad, temp3$dyad + idConvention)
	temp4 <- suppressMessages(plyr::join(datai, temp3))
	temp4$roleNew <- factor(temp4$role, levels=c("d0","d1"), labels=c(dist0name, dist1name)) 
			
	plotData <- temp4[stats::complete.cases(temp4), ]	
	plotTitle <- as.character(unique(datai$dyad))
						
	plots[[i]] <- ggplot(plotData, aes_string(x="time")) +
	  geom_line(aes_string(y= "obs_deTrend", color="roleNew"), linetype="dotted", size= .8, na.rm=T) +
	  geom_line(aes_string(y="pred", color="roleNew"), size= .8, na.rm=T) + 
	  scale_color_manual(name="Role", values=c("red","blue")) +
	  ylab(plot_obs_name) +
	  ylim(min, max) +
	  annotate("text", x=-Inf, y=-Inf, hjust=0, vjust=0, label="Dots = Observed; Lines = Predicted", size=3) +
	  labs(title= "Dyad ID:", subtitle= plotTitle) +
	  theme(plot.title=element_text(size=11)) +
	  theme(plot.subtitle=element_text(size=10))			
  }
  
  if(printPlots==T){print(plots)}	
  return(plots)
}

###################### cloResids

#' Produces histograms of the residuals from the oscillator model for each dyad.
#' 
#' @param derivData A dataframe that was produced with the "estDerivs" function.
#' @param whichModel Whether the model to be estimated is the uncoupled-oscillator ("uncoupled") or the coupled-oscillator ("coupled").
#' @param printPlots If true (the default) plots are displayed on the screen.
#' @examples
#' # See vignettes for examples.
#' 
#' @return A list with the histograms of the residuals for each dyad.

#' @import ggplot2
#' @export

cloResids <- function(derivData, whichModel, printPlots=T)
{
  basedata <- derivData
  
  if(whichModel != "uncoupled" & whichModel != "coupled") {
  	stop("the model type must be either uncoupled or coupled")
	
	} else if (whichModel == "uncoupled"){
	  model <- stats::formula(d2 ~ dist0:obs_deTrend + dist0:d1 + dist1:obs_deTrend + dist1:d1 -1)

      } else if (whichModel == "coupled"){
      	model <- stats::formula(d2 ~ dist0:obs_deTrend + dist0:d1 + dist0:p_obs_deTrend + dist0:p_d1 + dist1:obs_deTrend + dist1:d1 + dist1:p_obs_deTrend + dist1:p_d1 -1)
  }	   

  newDiD <- unique(factor(basedata$dyad))
  plots <- list()
  resid <- list()
	
  for (i in 1:length(newDiD)){
	datai <- basedata[basedata$dyad == newDiD[i], ]
	m <- stats::lm(model, na.action=na.exclude, data=datai) 
	plotTitle <- as.character(unique(datai$dyad))
	resid[[i]] <- m$residuals
	plotResid <- data.frame(resid[[i]])
	colnames(plotResid) <- "Residuals"
						
	plots[[i]] <- ggplot(plotResid, aes_string(x="Residuals")) +
	geom_histogram(color="black", fill="grey") +
	labs(title= "Dyad ID:", subtitle= plotTitle) +
	theme(plot.title=element_text(size=11)) +
	theme(plot.subtitle=element_text(size=10))		
  }
	
  if(printPlots==T){print(plots)}
  return(plots)
}


#################### cloPlotTraj

#' Plots the bivariate state variable's clo model-predicted temporal trajectories for each latent profile of clo parameters.
#' 
#' @param prepData A dataframe that was produced with the "dataPrep" function.
#' @param paramEst A dataframe created by indivClo containing the clo parameter estimates for each dyad.
#' @param n_profiles The number of latent profiles.
#' @param dist0name An optional name for the level-0 of the distinguishing variable (e.g., "Women"). Default is dist0.
#' @param dist1name An optional name for the level-1 of the distinguishing variable (e.g., "Men"). Default is dist1
#' @param plot_obs_name An optional name for the observed state variable to appear on plots (e.g., "Emotional Experience").
#' @param minMax An optional vector with desired minimum and maximum quantiles to be used for setting the y-axis range on the plots, e.g., minMax <- c(.1, .9) would set the y-axis limits to the 10th and 90th percentiles of the observed state variables. If not provided, the default is to use the minimum and maximum observed values of the state variables.
#' @param time_length An optional value specifying how many time points to plot across. Default is the 75th percentile for the observed time variable.
#' @param printPlots If true (the default) plots are displayed on the screen.
#' @examples
#' # See vignettes for examples.
#' 
#' @return The function returns the plots as a list. 

#' @import ggplot2
#' @export

cloPlotTraj <- function(prepData, paramEst, n_profiles, dist0name=NULL, dist1name=NULL, plot_obs_name = NULL, minMax=NULL, time_length=NULL, printPlots=T)
{    
  if(is.null(dist0name)){dist0name <- "dist0"}
  if(is.null(dist1name)){dist1name <- "dist1"}
  if(is.null(plot_obs_name)){plot_obs_name <- "observed"}
  
  if(is.null(minMax)){
    min <- min(prepData$obs_deTrend, na.rm=T)
    max <- max(prepData$obs_deTrend, na.rm=T)
  } else {
    min <- stats::quantile(prepData$obs_deTrend, minMax[1], na.rm=T)
    max <- stats::quantile(prepData$obs_deTrend, minMax[2],  na.rm=T)
  }
  
  if(is.null(time_length)){time_length <- as.numeric(stats::quantile(prepData$time, prob=.75))}
  
  prepData <- prepData[stats::complete.cases(prepData), ]	
  paramEst <- paramEst[stats::complete.cases(paramEst), ]	
  
  vars1 <- c("obs_0","d1_0","p_obs_0","p_d1_0","obs_1","d1_1","p_obs_1","p_d1_1")
  temp1 <- paramEst[vars1]
  lpa <- mclust::Mclust(temp1, G=n_profiles)
  profileParams <- as.data.frame(lpa$parameters$mean) 
  
  plots <- list()
  
  for(i in 1:n_profiles){
    statedata0 <- prepData[prepData$dist0 == 1 & prepData$time ==1,] 
    start0 <- stats::median(statedata0$obs_deTrend, na.rm=T)
    statedata1 <- prepData[prepData$dist0 == 0 & prepData$time ==1,] 
    start1 <- stats::median(statedata1$obs_deTrend, na.rm=T)
    
    plotTimes <- seq(1, time_length, by=1)
    
    state <- c("y1"=start0, "y2"=0, "y3"=start1, "y4"=0)
    
    temp1 <- profileParams[ ,i]
    names <- rownames(profileParams)
    names(temp1) <- names
    paramsi <- temp1
    
    temp2 <- as.data.frame(deSolve::ode(y=state, times=plotTimes, func=cloCoupledOde, parms= paramsi))
    vars2 <- c("y2", "y4")
    temp3 <- temp2[ ,!(names(temp2) %in% vars2)]
    names(temp3) <- c("time","d0pred","d1pred")
    temp4 <- stats::reshape(temp3, direction='long', varying=c("d0pred","d1pred"), timevar="role", times=c("d0","d1"), v.names=c("pred"), idvar="time")
    temp4$roleNew <- factor(temp4$role, levels=c("d0","d1"), labels=c(dist0name, dist1name)) 
    
    plotData <- temp4[stats::complete.cases(temp4), ]	
    profileName <- paste("Profile", i , sep="_")
    
    plotsi <- ggplot(plotData, aes_string(x="time")) +
      geom_line(aes_string(y="pred", color="roleNew"), linetype="solid", size=1, na.rm=T) +
      scale_color_manual(name="Role", values=c("black","gray47")) +
      ylab(plot_obs_name) +
      ylim(min, max) +
      labs(title=profileName, subtitle= "Predicted Trajectory") +
      theme(plot.title=element_text(size=11))
    
    plots[[i]] <- plotsi
  }
  if(printPlots==T){print(plots)}
  return(plots) 
} 
