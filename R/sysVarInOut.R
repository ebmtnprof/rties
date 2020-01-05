
################### sysVarOut

#' Provides results for predicting the system variable from the latent profiles of the dynamic parameters. 
#' 
#' The system variable can be either dyadic (sysVarType = "dyadic"), where both partners have the same score (e.g., relationship length) or individual (sysVarType = "indiv"), where the partners can have different scores (e.g., age). For dyadic system variables, the only predictor is profile membership and the model is a regular regression model since all variables are at the level of the dyad. If the system variable is individual then the model is a random-intercept dyadic model and 3 models are estimated: 1) the main effect of profile membership, 2) main effects of profile membership and the distinguishing variable, and 3) the interaction of profile membership and the distinguishing variable. If the system variable is not normally distributed, any of the generalized linear models supported by glm (for dyadic system variables) or glmmPQL (for individual system variables) are available by specifying the "family" distribution.
#' 
#' @param fullData A dataframe created by the "makeFullData" function.
#' @param sysVar_name The name of the variable in the dataframe that contains the system variable to be predicted by profile membership. 
#' @param sysVarType Whether the system variable is "dyadic", which means both partners have the same score, or "indiv" which means the partners can have different scores
#' @param dist0name An optional name for the level-0 of the distinguishing variable to appear as plot labels (e.g., "Women"). Default is dist0.
#' @param dist1name An optional name for the level-1 of the distinguishing variable to appear as plot labels (e.g., "Men"). Default is dist1
#' @param plot_sysVar_name An optional name for the system variable to appear as plot labels (e.g., "Satisfaction"). Default is sysVar.
#' @param minMax An optional vector with desired minimum and maximum quantiles to be used for setting the y-axis range on the plots, e.g., minMax <- c(.1, .9) would set the y-axis limits to the 10th and 90th percentiles of the observed state variables. If not provided, the default is to use the minimum and maximum observed values of the state variables.
#' @param family An optional argument specifying the error distribution and link function to be used in the model. Any of the "family" options supported by glm (for dyadic system variables) or glmmPQL (for individual system variables) are available. Default is gaussian.
#' @param printPlots Controls whether or not the plots are printed. Default is "true".
#' 
#' @return For normally distributed system variables, the function returns a list including the lm or lme objects containing the full results for each model (called "models"). Similarly, for non-normal system variables, the function returns a list of the glm or glmmPQL objects containing the full results for the models. By default, the function also displays histograms of the residuals and plots of the predicted values against observed values for each model, but these can be turned off by setting printPlots=F. 

#' @export

sysVarOut <- function(fullData, sysVar_name, sysVarType, dist0name=NULL, dist1name=NULL, plot_sysVar_name=NULL, minMax=NULL, family=NULL, printPlots=T)
{
  basedata <- fullData
  
  if(is.null(dist0name)){dist0name <- "dist0"}
  if(is.null(dist1name)){dist1name <- "dist1"}
  if(is.null(plot_sysVar_name)){plot_sysVar_name <- "System_Variable"}
  if(is.null(family)){family <- "gaussian"}
	
  colnames(basedata)[colnames(basedata)== sysVar_name] <- "sysVar" 	

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
	
  basedata$dist1 <- ifelse(basedata$dist0 == 1, 0, 1)
  basedata$dist <- factor(basedata$dist0, labels=c(dist1name, dist0name))
	
  if (sysVarType == "dyadic"){	
	basedata <- basedata[!duplicated(basedata$dyad), ]	
	
	if (family == "gaussian"){
	  profile <- lm(sysVar ~ profile, data= basedata, na.action=na.exclude)
	  profilePred <- fitted(profile)
	} else {
	  profile <- glm(sysVar ~ profile, data= basedata, na.action=na.exclude, family=family)
	  profilePred <- fitted(profile)
    }
      
    if(printPlots==T & family=="gaussian"){
	  ylabName <- paste(plot_sysVar_name, "predicted", sep="_")
	  xlabName <- paste(plot_sysVar_name, "observed", sep="_")

	  hist(residuals(profile))
	  plot(profilePred ~ basedata$sysVar, xlim=c(min, max), ylim=c(min, max), ylab=ylabName, xlab=xlabName, main="Profile Model")
	 }
  }

  if (sysVarType == "indiv"){
	if (family == "gaussian"){
	
	profile <- nlme::lme(sysVar ~ profile, random= ~ 1 | dyad, data= basedata, na.action=na.exclude, control=nlme::lmeControl(opt="optim"), method="ML")
	
	profilePlusDist <- nlme::lme(sysVar ~ profile + dist, random= ~ 1 | dyad, data= basedata, na.action=na.exclude, control=nlme::lmeControl(opt="optim"), method="ML")

	profileByDist <- nlme::lme(sysVar ~ profile * dist, random= ~ 1 | dyad, data= basedata, na.action=na.exclude, control=nlme::lmeControl(opt="optim"), method="ML")
    } else {
	
	  profile <- MASS::glmmPQL(sysVar ~ profile, random= ~ 1 | dyad, data= basedata, na.action=na.exclude, control=nlme::lmeControl(opt="optim"), family=family)
	
	  profilePlusDist <- MASS::glmmPQL(sysVar ~ profile + dist, random= ~ 1 | dyad, data= basedata, na.action=na.exclude, control=nlme::lmeControl(opt="optim"), family=family)

	  profileByDist <- MASS::glmmPQL(sysVar ~ profile * dist, random= ~ 1 | dyad, data= basedata, na.action=na.exclude, control=nlme::lmeControl(opt="optim"), family=family)
	  }
	
	profilePred <- fitted(profile)
	profilePlusDistPred <- fitted(profilePlusDist)
	profileByDistPred <- fitted(profileByDist)
  	
    if(printPlots==T & family=="gaussian"){
	  ylabName <- paste(plot_sysVar_name, "predicted", sep="_")
	  xlabName <- paste(plot_sysVar_name, "observed", sep="_")
	  
	  hist(residuals(profile))
	  plot(profilePred ~ basedata$sysVar, xlim=c(min, max), ylim=c(min, max), ylab=ylabName, xlab=xlabName, main="Profile Model")
	
	  hist(residuals(profilePlusDist))
	  plot(profilePlusDistPred ~ basedata$sysVar, xlim=c(min, max), ylim=c(min, max), ylab=ylabName, xlab=xlabName, main="Profile Plus Dist Model")

	  hist(residuals(profileByDist))
	  plot(profileByDistPred ~ basedata$sysVar, xlim=c(min, max), ylim=c(min, max), ylab=ylabName, xlab=xlabName, main="Profile By Dist Model")
	  
	  basedata$dist <- factor(basedata$dist0, labels=c(dist1name, dist0name))
	  interact <- ggplot(basedata, aes(x=profile, y=sysVar, fill=dist)) +
                    geom_boxplot() + 
                    scale_fill_manual(values=c("gray88","gray60")) + 
                    ylab(plot_sysVar_name)
      print(interact)
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
#' If there are 2 profiles, then binomial regression models are used. If there are more than 2 profiles then multinomial regression is used. The system variable can be either dyadic (sysVarType = "dyadic"), where both partners have the same score (e.g., relationship length) or individual (sysVarType = "indiv"), where the partners can have different scores (e.g., age). For dyadic system variables, a couple's shared score is the only predictor of their profile membership (called "sysVar"). For individual system variables, two models are tested, one with the main effects of both partner's system variable ("sysVarMain") and one with the main effects and their interaction ("sysVarInteract"). In both cases an intercept-only model is included as a comparison point (called "base"). The function returns a list of the full model results and produces plots of the probability of profile membership against the system variable(s).
#' 
#' @param fullData A dataframe created by the makeFullData function.
#' @param sysVar_name The name of the variable in the dataframe that contains the system variable to be predicted by profile membership. 
#' @param sysVarType Whether the system variable is "dyadic", which means both partners have the same score, or "indiv" which means the partners can have different scores
#' @param n_profiles The number of latent profiles.
#' @param dist0name An optional name for the level-0 of the distinguishing variable (e.g., "Women"). Default is dist0.
#' @param dist1name An optional name for the level-1 of the distinguishing variable (e.g., "Men"). Default is dist1
#' @param plot_sysVar_name An optional name for the system variable being predicted (e.g., "Satisfaction"). Default is sysVar.
#' 
#' @return A list of model results and, by default, plots of profile membership against the system variable(s), but these can be turned off by setting printPlots=F.

#' @export

sysVarIn <- function(fullData, sysVar_name, sysVarType, n_profiles, dist0name=NULL, dist1name=NULL, plot_sysVar_name=NULL){

  basedata <- fullData
  
  if(is.null(dist0name)){dist0name <- "dist0"}
  if(is.null(dist1name)){dist1name <- "dist1"}
  if(is.null(plot_sysVar_name)){plot_sysVar_name <- "sysVar"}
  
  if(sysVarType != "indiv" & sysVarType != "dyadic") {stop("the sysVarType must be either indiv or dyadic")}

  colnames(basedata)[colnames(basedata)== sysVar_name] <- "sysVar" 	

  if(sysVarType == "dyadic"){
    
    basedata <- basedata[!duplicated(basedata$dyad), ]
    
    if(n_profiles == 2){
      base <- glm(profileN ~ 1, data=basedata, na.action=na.exclude, family="binomial")
      sysVarMain <- glm(profileN ~ sysVar, data=basedata, na.action=na.exclude, family="binomial")    	
     } else {
     	 base <- nnet::multinom(profileN ~ 1, data=basedata, na.action=na.exclude,)
     	 sysVarMain <- nnet::multinom(profileN ~ sysVar, data=basedata, na.action=na.exclude)
     }	
	  	plot(basedata$sysVar, basedata$profile)
  }
  
    if(sysVarType == "indiv"){
    
    data1 <- subset(basedata, select=c(dyad, sysVar, dist0, profileN))
    data2 <-  stats::reshape(data1, idvar="dyad", timevar = "dist0", direction= "wide")   
    data3 <- dplyr::rename(data2, dyad=dyad, sysVar0=sysVar.1, profileN1= profileN.1, sysVar1= sysVar.0, profileN=profileN.0)
    data4 <- subset(data3, select=-c(profileN1))   
    basedata <- data4
    
    sysVar0name <- paste(plot_sysVar_name, dist0name, sep="_")
	sysVar1name <- paste(plot_sysVar_name, dist1name, sep="_")
	sysVar01name <- paste(sysVar0name, sysVar1name, sep=":")
   
    if(n_profiles == 2){
      base <- glm(profileN ~ 1, data=basedata, na.action=na.exclude, family="binomial")
      sysVarMain <- glm(profileN ~ sysVar0 + sysVar1, data=basedata, na.action=na.exclude, family="binomial") 
      names(sysVarMain$coefficients) <- c("Intercept", sysVar0name, sysVar1name)
      sysVarInteract <- glm(profileN ~ sysVar0 * sysVar1, data=basedata, na.action=na.exclude, family="binomial") 
      names(sysVarInteract$coefficients) <- c("Intercept", sysVar0name, sysVar1name, sysVar01name)
      
      if(is.factor(basedata$sysVar0)){
	  	print(interactions::cat_plot(sysVarInteract, pred=sysVar0, modx=sysVar1, y.label="Prob Profile = 2", x.label=sysVar0name, legend.main=sysVar1name, colors="Greys", interval=T))
	   } else {
	  	print(interactions::interact_plot(sysVarInteract, pred=sysVar0, modx=sysVar1, y.label="Prob Profile = 2", x.label=sysVar0name, legend.main=sysVar1name, colors="Greys", interval=T))
	   }

    } else {
    	base <- nnet::multinom(profileN ~ 1, data=basedata)
    	sysVarMain <- nnet::multinom(profileN ~ sysVar0 + sysVar1, data=basedata, na.action=na.exclude)
    	names(sysVarMain$coefnames) <- c("Intercept", sysVar0name, sysVar1name)
    	sysVarInteract <- nnet::multinom(profileN ~ sysVar0 * sysVar1, data=basedata, na.action=na.exclude)
    	names(sysVarInteract$coefnames) <- c("Intercept", sysVar0name, sysVar1name, sysVar01name)
 
    sysVar0L <- mean(basedata$sysVar0) - sd(basedata$sysVar0)
    sysVar0H <- mean(basedata$sysVar0) + sd(basedata$sysVar0)
    sysVar1L <- mean(basedata$sysVar1) - sd(basedata$sysVar1)
    sysVar1H <- mean(basedata$sysVar1) + sd(basedata$sysVar1)

    dataTemp<- matrix(c(sysVar0L, sysVar0H, sysVar0L, sysVar0H, sysVar1L, sysVar1L, sysVar1H, sysVar1H), nrow=4, ncol=2)
    dataTemp2 <- data.frame(dataTemp)
    colnames(dataTemp2) <- c("sysVar0", "sysVar1")

    prob <- data.frame(predict(sysVarInteract, newdata=dataTemp2, "probs"))
    prob$sysVar0 <- c(1,2,1,2)
    prob$sysVar1 <- c(1,1,2,2)
    prob$sysVar0 <- factor(prob$sysVar0, levels=c(1,2), labels=c("Low", "High"))
    prob$sysVar1 <- factor(prob$sysVar1, levels=c(1,2), labels=c("Low", "High"))

    if(n_profiles==3){
      prob1 <- subset(prob, select=c(X0, sysVar0, sysVar1))
      prob2 <- subset(prob, select=c(X1, sysVar0, sysVar1))
      prob3 <- subset(prob, select=c(X2, sysVar0, sysVar1))

      p1 <- ggplot(prob1, aes(x = sysVar0, y=X0, group=sysVar1)) + 
        geom_line(aes(linetype=sysVar1)) +
        ylim(0,1) +
        labs(title="Profile-1",y="Probabilty", x=sysVar0name) + 
        scale_linetype_discrete(name=sysVar1name)
  
      p2 <- ggplot(prob2, aes(x = sysVar0, y=X1, group=sysVar1)) + 
        geom_line(aes(linetype=sysVar1)) +
        ylim(0,1) +
        labs(title="Profile-2",y="Probabilty", x=sysVar0name) + 
        scale_linetype_discrete(name=sysVar1name)
  
      p3 <- ggplot(prob3, aes(x = sysVar0, y=X2, group=sysVar1)) + 
        geom_line(aes(linetype=sysVar1)) +
        ylim(0,1) +
        labs(title="Profile-3",y="Probabilty", x=sysVar0name) + 
        scale_linetype_discrete(name=sysVar1name)

      pAll <- gridExtra::grid.arrange(p1,p2,p3, nrow=3)
     }

     if(n_profiles==4){
       prob1 <- subset(prob, select=c(X0, sysVar0, sysVar1))
       prob2 <- subset(prob, select=c(X1, sysVar0, sysVar1))
       prob3 <- subset(prob, select=c(X2, sysVar0, sysVar1))
       prob4 <- subset(prob, select=c(X3, sysVar0, sysVar1))

       p1 <- ggplot(prob1, aes(x = sysVar0, y=X0, group=sysVar1)) + 
         geom_line(aes(linetype=sysVar1)) +
         ylim(0,1) +
         labs(title="Profile-1",y="Probabilty", x=sysVar0name) + 
         scale_linetype_discrete(name=sysVar1name)
  
       p2 <- ggplot(prob2, aes(x = sysVar0, y=X1, group=sysVar1)) + 
         geom_line(aes(linetype=sysVar1)) +
         ylim(0,1) +
         labs(title="Profile-2",y="Probabilty", x=sysVar0name) + 
         scale_linetype_discrete(name=sysVar1name)
  
       p3 <- ggplot(prob3, aes(x = sysVar0, y=X2, group=sysVar1)) + 
         geom_line(aes(linetype=sysVar1)) +
         ylim(0,1) +
         labs(title="Profile-3",y="Probabilty", x=sysVar0name) + 
         scale_linetype_discrete(name=sysVar1name)
  
       p4 <- ggplot(prob4, aes(x = sysVar0, y=X3, group=sysVar1)) + 
         geom_line(aes(linetype=sysVar1)) +
         ylim(0,1) +
         labs(title="Profile-4",y="Probabilty", x=sysVar0name) + 
         scale_linetype_discrete(name=sysVar1name)

         pAll <- gridExtra::grid.arrange(p1,p2,p3,p4, nrow=4)
       }

       if(n_profiles > 4) {message("plots are not provided if there are more than 4 profiles")}

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

