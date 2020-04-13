
################### sysVarOut

#' Provides results for predicting the system variable from the latent profiles of the dynamic parameters. 
#' 
#' The system variable can be either dyadic (sysVarType = "dyadic"), where both partners have the same score (e.g., relationship length) or individual (sysVarType = "indiv"), where the partners can have different scores (e.g., age). For dyadic system variables, the only predictor is profile membership and the model is a regular regression model since all variables are at the level of the dyad. If the system variable is individual then the model is a random-intercept dyadic model and 3 models are estimated: 1) the main effect of profile membership, 2) main effects of profile membership and the distinguishing variable, and 3) the interaction of profile membership and the distinguishing variable. If the system variable is not normally distributed, any of the generalized linear models supported by glm (for dyadic system variables) or glmmPQL (for individual system variables) are available by specifying the "family" distribution.
#' 
#' @param fullData A dataframe created by the "makeFullData" function.
#' @param sysVar_name The name of the variable in the dataframe that contains the system variable to be predicted by profile membership. 
#' @param sysVarType Whether the system variable is "dyadic", which means both partners have the same score, or "indiv" which means the partners can have different scores
#' @param dist0name An optional name for the level-0 of the distinguishing variable (e.g., "Women"). Default is dist0.
#' @param dist1name An optional name for the level-1 of the distinguishing variable (e.g., "Men"). Default is dist1
#' @param family An optional argument specifying the error distribution and link function to be used in the model. Any of the "family" options supported by glm (for dyadic system variables) or glmmPQL (for individual system variables) are available. Default is gaussian.
#' 
#' @return For normally distributed system variables, the function returns a list including the lm or lme objects containing the full results for each model (called "models"). Similarly, for non-normal system variables, the function returns a list of the glm or glmmPQL objects containing the full results for the models.  

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
	
  if (sysVarType == "dyadic"){	
	  basedata <- basedata[!duplicated(basedata$dyad), ]	
	
	  if (family == "gaussian"){
	    profile <- stats::lm(sysVar ~ profile, data= basedata, na.action=na.exclude)
	    profilePred <- stats::fitted(profile)
	  } else {
	    profile <- stats::glm(sysVar ~ profile, data= basedata, na.action=na.exclude, family=family)
	    profilePred <- stats::fitted(profile)
    }
  }

  if (sysVarType == "indiv"){
	  
    if (family == "gaussian"){
	    profile <- nlme::lme(sysVar ~ profile, random= ~ 1 | dyad, data= basedata, na.action=na.exclude, control=nlme::lmeControl(opt="optim"), method="ML")
	    profilePlusDist <- nlme::lme(sysVar ~ profile + dist, random= ~ 1 | dyad, data= basedata, na.action=na.exclude, control=nlme::lmeControl(opt="optim"), method="ML")
	    profileByDist <- nlme::lme(sysVar ~ profile * dist, random= ~ 1 | dyad, data= basedata, na.action=na.exclude, control=nlme::lmeControl(opt="optim"), method="ML")
    } else {
	    profile <- MASS::glmmPQL(sysVar ~ profile, random= ~ 1 | dyad, data= basedata, na.action=na.exclude, control=nlme::lmeControl(opt="optim"), family=family)=
	    profilePlusDist <- MASS::glmmPQL(sysVar ~ profile + dist, random= ~ 1 | dyad, data= basedata, na.action=na.exclude, control=nlme::lmeControl(opt="optim"), family=family)
	    profileByDist <- MASS::glmmPQL(sysVar ~ profile * dist, random= ~ 1 | dyad, data= basedata, na.action=na.exclude, control=nlme::lmeControl(opt="optim"), family=family)
	  }
	}

  if(sysVarType == "dyadic"){
	  models <- list(profile=profile)
	  message("Model name is profile")
  }
	
  if(sysVarType == "indiv"){
	  models <- list(profile=profile, profilePlusDist=profilePlusDist, profileByDist=profileByDist)
	  message("Model names are profile, profilePlusDist and profileByDist")
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
#' @param dist0name An optional name for the level-0 of the distinguishing variable (e.g., "Women"). Default is dist0.
#' @param dist1name An optional name for the level-1 of the distinguishing variable (e.g., "Men"). Default is dist1
#' 
#' @return A list including the glm or multinom objects containing the full results for each model (called "models"). 

#' @export

sysVarIn <- function(fullData, sysVar_name, sysVarType, n_profiles, dist0name=NULL, dist1name=NULL){

  basedata <- fullData
  
  if(is.null(dist0name)){dist0name <- "dist0"}
  if(is.null(dist1name)){dist1name <- "dist1"}
  
  sysVar0name <- paste(sysVar_name, dist0name, sep="_")
  sysVar1name <- paste(sysVar_name, dist1name, sep="_")
  sysVar01name <- paste(sysVar0name, sysVar1name, sep=":")
  
  if(sysVarType != "indiv" & sysVarType != "dyadic") {stop("the sysVarType must be either indiv or dyadic")}

  colnames(basedata)[colnames(basedata)== sysVar_name] <- "sysVar" 	

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
      basedata <- data3
    
        if(n_profiles == 2){
          base <- stats::glm(profileN ~ 1, data=basedata, na.action=na.exclude, family="binomial")
          sysVarMain <- stats::glm(profileN ~ sysVar0 + sysVar1, data=basedata, na.action=na.exclude, family="binomial") 
          names(sysVarMain$coefficients) <- c("Intercept", sysVar0name, sysVar1name)
          sysVarInteract <- stats::glm(profileN ~ sysVar0 * sysVar1, data=basedata, na.action=na.exclude, family="binomial") 
          names(sysVarInteract$coefficients) <- c("Intercept", sysVar0name, sysVar1name, sysVar01name)
        } else {
    	     base <- nnet::multinom(profileN ~ 1, data=basedata)
    	     sysVarMain <- nnet::multinom(profileN ~ sysVar0 + sysVar1, data=basedata, na.action=na.exclude)
    	     names(sysVarMain$coefnames) <- c("Intercept", sysVar0name, sysVar1name)
    	     sysVarInteract <- nnet::multinom(profileN ~ sysVar0 * sysVar1, data=basedata, na.action=na.exclude)
    	     names(sysVarInteract$coefnames) <- c("Intercept", sysVar0name, sysVar1name, sysVar01name)
        }
    }
  
    if(sysVarType == "dyadic"){
	    models <- list(base=base, sysVarMain=sysVarMain)
	    message("Model names are base and sysVarMain")
    }
	
    if(sysVarType == "indiv"){
	    models <- list(base=base, sysVarMain=sysVarMain, sysVarInteract=sysVarInteract)
	    ("Model names are base, sysVarMain and sysVarInteract")
    }
	
  output <- list(models=models)
}

