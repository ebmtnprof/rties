
################### sysVarOut

#' Provides results for predicting the system variable from the latent profiles of the dynamic parameters. 
#' 
#' The system variable can be either dyadic (sysVarType = "dyadic"), where both partners have the same score (e.g., relationship length) or individual (sysVarType = "indiv"), where the partners can have different scores (e.g., age). For dyadic system variables, the only predictor is profile membership and the model is a regular regression model since all variables are at the level of the dyad. If the system variable is individual then the model is a random-intercept dyadic model and 3 models are estimated: 1) the main effect of profile membership, 2) main effects of profile membership and the distinguishing variable, and 3) the interaction of profile membership and the distinguishing variable. If the system variable is not normally distributed, any of the generalized linear models supported by glm (for dyadic system variables) or glmer (for individual system variables) are available by specifying the "family" distribution.
#' 
#' @param fullData A dataframe created by the "makeFullData" function.
#' @param sysVar_name The name of the variable in the dataframe that contains the system variable to be predicted by profile membership. 
#' @param sysVarType Whether the system variable is "dyadic", which means both partners have the same score, or "indiv" which means the partners can have different scores
#' @param dist0name An optional name for the level-0 of the distinguishing variable (e.g., "Women"). Default is dist0.
#' @param dist1name An optional name for the level-1 of the distinguishing variable (e.g., "Men"). Default is dist1
#' @param family An optional argument specifying the error distribution and link function to be used in the model. Any of the "family" options supported by glm (for dyadic system variables) or glmer (for individual system variables) are available. Default is gaussian.
#' @examples
#' data <- rties_ExampleDataShort
#' newData <- dataPrep(basedata=data, dyadId="couple", personId="person", 
#' obs_name="dial", dist_name="female", time_name="time", time_lag=2)
#' ic <- indivInertCoord(prepData=newData, whichModel="inertCoord")
#' profiles <- inspectProfiles(whichModel="inertCoord", prepData=newData, 
#' paramEst=ic$params, n_profiles=2)
#' fullData <- makeFullData(basedata=data, dyadId="couple", personId="person", 
#' dist_name="female", lpaData=profiles, params=ic$params)
#' sysOut <- sysVarOut(fullData=fullData, sysVar_name="conflict", sysVarType="indiv")
#' summary(sysOut$models$profilePlusDist)
#' 
#' @return For normally distributed system variables, the function returns a list including the lm or lme objects containing the full results for each model (called "models"). Similarly, for non-normal system variables, the function returns a list of the glm or glmer objects containing the full results for the models.  

#' @export

sysVarOut <- function(fullData, sysVar_name, sysVarType, dist0name=NULL, dist1name=NULL, family=NULL)
{
  basedata <- fullData 

  if(is.null(family)){family <- "gaussian"}
  if(is.null(dist0name)){dist0name <- "dist0"}
  if(is.null(dist1name)){dist1name <- "dist1"}
  
  colnames(basedata)[colnames(basedata)== sysVar_name] <- "sysVar" 	

  if(sysVarType != "indiv" & sysVarType != "dyadic") {
	stop("the sysVarType must be either indiv or dyadic")
  }
	
  basedata$dist1 <- ifelse(basedata$dist0 == 1, 0, 1)
  basedata$dist <- factor(basedata$dist0, labels=c(dist1name, dist0name))
  basedata <- basedata[stats::complete.cases(basedata), ]
	
  if (sysVarType == "dyadic"){	
	  basedata <- basedata[!duplicated(basedata$dyad), ]
	  basedata <- basedata[stats::complete.cases(basedata), ]
	
	  if (family == "gaussian"){
	    base <- stats::lm(sysVar ~ 1, data= basedata, na.action=na.exclude)
	    profile <- stats::lm(sysVar ~ profile, data= basedata, na.action=na.exclude)
	  } else {
	    base <- stats::glm(sysVar ~ 1, data= basedata, na.action=na.exclude, family=family)
	    profile <- stats::glm(sysVar ~ profile, data= basedata, na.action=na.exclude, family=family)
    }
  }

  if (sysVarType == "indiv"){
	  
    if (family == "gaussian"){
	    base <- nlme::lme(sysVar ~ 1, random= ~ 1 | dyad, data= basedata, na.action=na.exclude, control=nlme::lmeControl(opt="optim"), method="ML")
      profile <- nlme::lme(sysVar ~ profile, random= ~ 1 | dyad, data= basedata, na.action=na.exclude, control=nlme::lmeControl(opt="optim"), method="ML")
	    profilePlusDist <- nlme::lme(sysVar ~ profile + dist, random= ~ 1 | dyad, data= basedata, na.action=na.exclude, control=nlme::lmeControl(opt="optim"), method="ML")
	    profileByDist <- nlme::lme(sysVar ~ profile * dist, random= ~ 1 | dyad, data= basedata, na.action=na.exclude, control=nlme::lmeControl(opt="optim"), method="ML")
    } else {
      base <- lme4::glmer(sysVar ~ 1 + (1 | dyad), data= basedata, na.action=na.exclude, control=lme4::glmerControl(), family=family)
      profile <- lme4::glmer(sysVar ~ profile + (1 | dyad), data= basedata, na.action=na.exclude, control=lme4::glmerControl(), family=family)
	    profilePlusDist <- lme4::glmer(sysVar ~ profile + dist + (1 | dyad), data= basedata, na.action=na.exclude, control=lme4::glmerControl(), family=family)
	    profileByDist <- lme4::glmer(sysVar ~ profile * dist + (1 | dyad), data= basedata, na.action=na.exclude, control=lme4::glmerControl(), family=family)
	  }
	}

  if(sysVarType == "dyadic"){
	  models <- list(base=base, profile=profile)
	  message("Model names are base & profile")
  }
	
  if(sysVarType == "indiv"){
	  models <- list(base=base, profile=profile, profilePlusDist=profilePlusDist, profileByDist=profileByDist)
	  message("Model names are base, profile, profilePlusDist and profileByDist")
  }
	output <- list(models=models)
}

################### sysVarIn

#' Provides results for predicting couples' latent profile membership from the system variable. 
#' 
#' If there are 2 profiles, then glm binomial regression models are used. If there are more than 2 profiles then multinomial regression is used (from the nnet package). The system variable can be either dyadic (sysVarType = "dyadic"), where both partners have the same score (e.g., relationship length) or individual (sysVarType = "indiv"), where the partners can have different scores (e.g., age). For dyadic system variables, a couple's shared score is the only predictor of their profile membership (called "sysVar"). For individual system variables, two models are tested, one with the main effects of both partner's system variable ("sysVarMain") and one with the main effects and their interaction ("sysVarInteract"). In both cases an intercept-only model is included as a comparison point (called "base"). The function returns a list of the full model results.
#' 
#' @param fullData A dataframe created by the makeFullData function.
#' @param sysVar_name The name of the variable in the dataframe that contains the system variable to be predicted by profile membership.
#' @param sysVarType Whether the system variable is "dyadic", which means both partners have the same score, or "indiv" which means the partners can have different scores
#' @param n_profiles The number of latent profiles.
#' @examples
#' data <- rties_ExampleDataShort
#' newData <- dataPrep(basedata=data, dyadId="couple", personId="person", 
#' obs_name="dial", dist_name="female", time_name="time", time_lag=2)
#' ic <- indivInertCoord(prepData=newData, whichModel="inertCoord")
#' profiles <- inspectProfiles(whichModel="inertCoord", prepData=newData, 
#' paramEst=ic$params, n_profiles=2)
#' fullData <- makeFullData(basedata=data, dyadId="couple", personId="person", 
#' dist_name="female", lpaData=profiles, params=ic$params)
#' sysIn <- sysVarIn(fullData=fullData, sysVar_name="conflict", sysVarType="indiv", n_profiles=2)
#' summary(sysIn$models$sysVarMain)
#' 
#' @return A list including the glm or multinom objects containing the full results for each model (called "models"). 

#' @export

sysVarIn <- function(fullData, sysVar_name, sysVarType, n_profiles){

  basedata <- fullData
  
  if(sysVarType != "indiv" & sysVarType != "dyadic") {stop("the sysVarType must be either indiv or dyadic")}

  colnames(basedata)[colnames(basedata)== sysVar_name] <- "sysVar" 	
  basedata <- basedata[stats::complete.cases(basedata), ]

  if(sysVarType == "dyadic"){
    
    basedata <- basedata[!duplicated(basedata$dyad), ]
    
    if(n_profiles == 2){
      base <- stats::glm(profileN ~ 1, data=basedata, na.action=na.exclude, family="binomial")
      sysVarMain <- stats::glm(profileN ~ sysVar, data=basedata, na.action=na.exclude, family="binomial")    	
     } else {
     	 base <- nnet::multinom(profileN ~ 1, data=basedata, na.action=na.exclude,)
     	 sysVarMain <- nnet::multinom(profileN ~ sysVar, data=basedata, na.action=na.exclude)
     }	
  }
  
    if(sysVarType == "indiv"){
      vars1 <- c("dyad", "sysVar", "dist0", "profileN")
      data1 <- basedata[vars1]
      data2 <-  stats::reshape(data1, idvar="dyad", timevar = "dist0", direction= "wide")   
      dyad <- sysVar.1 <- profileN.1 <- sysVar.0 <- profileN.0 <- NULL
      data3 <- dplyr::rename(data2, dyad=dyad, sysVar0=sysVar.1, profileN1= profileN.1, sysVar1= sysVar.0, profileN=profileN.0)
      basedata <- data3[stats::complete.cases(data3), ]
    
        if(n_profiles == 2){
          base <- stats::glm(profileN ~ 1, data=basedata, na.action=na.exclude, family="binomial")
          sysVarMain <- stats::glm(profileN ~ sysVar0 + sysVar1, data=basedata, na.action=na.exclude, family="binomial") 
          sysVarInteract <- stats::glm(profileN ~ sysVar0 * sysVar1, data=basedata, na.action=na.exclude, family="binomial") 
        } else {
    	     base <- nnet::multinom(profileN ~ 1, data=basedata)
    	     sysVarMain <- nnet::multinom(profileN ~ sysVar0 + sysVar1, data=basedata, na.action=na.exclude)
    	     sysVarInteract <- nnet::multinom(profileN ~ sysVar0 * sysVar1, data=basedata, na.action=na.exclude)
        }
    }
  
    if(sysVarType == "dyadic"){
	    models <- list(base=base, sysVarMain=sysVarMain)
	    message("Model names are base and sysVarMain")
    }
	
    if(sysVarType == "indiv"){
	    models <- list(base=base, sysVarMain=sysVarMain, sysVarInteract=sysVarInteract)
	    message("Model names are base, sysVarMain and sysVarInteract")
    }
	
  output <- list(models=models)
}


############### sysVarInResults

#' Produces results from sysVarIn.
#'
#' @param baseModel The name of the model that was produced by sysVarIn to be used as the null model for comparison (e.g., sysIn$models$base).
#' @param testModel The name of the model that was produced by sysVarIn that you want results for (e.g., sysIn$models$sysVarMain or sysIn$models$sysVarInteract).
#' @param n_profiles The number of latent profiles.
#' #' @examples
#' data <- rties_ExampleDataShort
#' newData <- dataPrep(basedata=data, dyadId="couple", personId="person", 
#' obs_name="dial", dist_name="female", time_name="time", time_lag=2)
#' ic <- indivInertCoord(prepData=newData, whichModel="inertCoord")
#' profiles <- inspectProfiles(whichModel="inertCoord", prepData=newData, 
#' paramEst=ic$params, n_profiles=2)
#' fullData <- makeFullData(basedata=data, dyadId="couple", personId="person", 
#' dist_name="female", lpaData=profiles, params=ic$params)
#' sysIn <- sysVarIn(fullData=fullData, sysVar_name="conflict", sysVarType="indiv", n_profiles=2)
#' sysVarInResults(baseModel=sysIn$models$base, testModel=sysIn$models$sysVarMain, n_profiles=2)
#' 
#' @return A list of results including a comparison of the test model to the null (either a LRT or Chisquare test depending on the model), a summary of the parameter estimates, exponentiated parameter estimates (e.g., odds ratios), and p values for the parameter estimates.

#' @export

sysVarInResults <- function(baseModel, testModel, n_profiles){
  
  if(n_profiles == 2){
    modelCompare <- stats::anova(baseModel, testModel, test="LRT")
    paramEst <- summary(testModel)
    oddsRatio <- exp(stats::coef(testModel))
    results <- list(modelCompare=modelCompare, paramEst=paramEst, oddsRatio=oddsRatio)
  }
  
  if(n_profiles > 2){
    modelCompare <- stats::anova(baseModel, testModel, test="Chisq")
    paramEst <- summary(testModel)
    z <- summary(testModel)$coefficients/summary(testModel)$standard.errors
    p <- (1 - stats::pnorm(abs(z), 0, 1)) * 2
    p <- round(p, 3)
    oddsRatio <- exp(stats::coef(testModel))
    results <- list(modelCompare=modelCompare, paramEst=paramEst, p=p, oddsRatio=oddsRatio)
  }
  return(results)
}

############### sysVarOutResults

#' Produces results from sysVarOut.
#'
#' @param baseModel The name of the model that was produced by sysVarOut to be used as the null model for comparison (e.g., sysOut$models$base).
#' @param testModel The name of the model that was produced by sysVarOut that you want results for (e.g., sysOut$models$profile, sysOut$models$profilePlusDist, sysOut$models$profileByDist).
#' @param Gaussian Whether the system variable is Gaussian. Default is true.
#' @examples
#' data <- rties_ExampleDataShort
#' newData <- dataPrep(basedata=data, dyadId="couple", personId="person", 
#' obs_name="dial", dist_name="female", time_name="time", time_lag=2)
#' ic <- indivInertCoord(prepData=newData, whichModel="inertCoord")
#' profiles <- inspectProfiles(whichModel="inertCoord", prepData=newData, 
#' paramEst=ic$params, n_profiles=2)
#' fullData <- makeFullData(basedata=data, dyadId="couple", personId="person", 
#' dist_name="female", lpaData=profiles, params=ic$params)
#' sysOut <- sysVarOut(fullData=fullData, sysVar_name="conflict", sysVarType="indiv")
#' sysVarOutResults(baseModel=sysOut$models$base, testModel=sysOut$models$profileByDist)
#' 
#' @return A list of results including an LRT comparison of the test model to the null, an omnibus anova test for the parameters in the model (this is identical to the LRT test for Gaussian dyadic system variables), a summary of the parameter estimates, and exponentiated parameter estimates (e.g., odds ratios) if Gaussian = FALSE.

#' @export
sysVarOutResults <- function(baseModel, testModel, Gaussian=TRUE){
  
  if(Gaussian == T){
    modelCompare <- stats::anova(baseModel, testModel)
    omnibus <- stats::anova(testModel)
    paramEst <- summary(testModel)
    results <- list(modelCompare=modelCompare, omnibus=omnibus, paramEst=paramEst)
  } else {
    modelCompare <- stats::anova(baseModel, testModel, test="LRT")
    paramEst <- summary(testModel)
    oddsRatio <- exp(summary(testModel)$coefficients[ ,1])
    results <- list(modelCompare=modelCompare, paramEst=paramEst, oddsRatio=oddsRatio)
  }
  return(results)
}
