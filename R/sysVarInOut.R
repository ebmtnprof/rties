
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
  basedata <- basedata[complete.cases(basedata), ]
	
  if (sysVarType == "dyadic"){	
	  basedata <- basedata[!duplicated(basedata$dyad), ]	
	
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
#' 
#' @return A list including the glm or multinom objects containing the full results for each model (called "models"). 

#' @export

sysVarIn <- function(fullData, sysVar_name, sysVarType, n_profiles){

  basedata <- fullData
  
  if(sysVarType != "indiv" & sysVarType != "dyadic") {stop("the sysVarType must be either indiv or dyadic")}

  colnames(basedata)[colnames(basedata)== sysVar_name] <- "sysVar" 	
  basedata <- basedata[complete.cases(basedata), ]

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
      basedata <- data3[complete.cases(data3), ]
    
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


##############################################

############### sysVarInPlots

#' Produces plots for interpreting the results from sysVarIn.
#'
#' @param fullData A dataframe created by the "makeFullData" function.
#' @param sysVar_name The name of the variable in the dataframe that contains the system variable.
#' @param sysVarType Whether the system variable is "dyadic", which means both partners have the same score, or "indiv" which means the partners can have different scores
#' @param n_profiles The number of latent profiles.
#' @param model The name of the model that is being interpreted (e.g., sysIn$models$sysVarInteract). Only needed when the system variable is "indiv" (e.g., individual scores for each partner)
#' @param dist0name An optional name for the level-0 of the distinguishing variable (e.g., "Women"). Default is dist0.
#' @param dist1name An optional name for the level-1 of the distinguishing variable (e.g., "Men"). Default is dist1
#' @param printPlots If true (the default) plots are displayed on the screen.
#' 
#' @return Single plots or a list of plots (depending on the model that is being interpreted).

#' @import ggplot2
#' @export

sysVarInPlots <- function(fullData, sysVar_name, sysVarType, n_profiles, model=NULL, dist0name=NULL, dist1name=NULL, printPlots=T){
  
  basedata <- fullData
  basedata <- basedata[complete.cases(basedata), ]
  colnames(basedata)[colnames(basedata)== sysVar_name] <- "sysVar" 	
  if(n_profiles > 4) {message("plots are not provided if there are more than 4 profiles") }
  
  if(is.null(dist0name)){dist0name <- "dist0"}
  if(is.null(dist1name)){dist1name <- "dist1"}
  
  if(sysVarType == "dyadic"){
    ### Dyadic 
    basedata <- basedata[!duplicated(basedata$dyad), ]	
    
    if(is.factor(basedata$sysVar)){
      pAll <- ggplot(basedata) +
        geom_bar(aes(x = sysVar, y = ..prop.., group = profile)) +
        facet_wrap(~ profile) +
        xlab(sysVar_name) +
        ylab("Proportion in Each Profile")
    } else {
      pAll <- ggplot(basedata, aes(x=sysVar, y=profileN)) +
        geom_point() +
        xlab(sysVar_name) +
        ylab("Profile")
    }
  }
  
  if(sysVarType == "indiv"){
    ### Indiv
    vars1 <- c("dyad", "sysVar", "dist0", "profileN")
    data1 <- basedata[vars1]
    data2 <-  stats::reshape(data1, idvar="dyad", timevar = "dist0", direction= "wide")   
    dyad <- sysVar.1 <- profileN.1 <- sysVar.0 <- profileN.0 <- NULL
    data3 <- dplyr::rename(data2, dyad=dyad, sysVar0=sysVar.1, profileN1= profileN.1, sysVar1= sysVar.0, profileN=profileN.0)
    basedata <- data3[complete.cases(data3), ]
    
    sysVar0name <- paste(sysVar_name, dist0name, sep="_")
    sysVar1name <- paste(sysVar_name, dist1name, sep="_")
    
    if(n_profiles == 2){
      ### 2 profiles
      if(is.factor(basedata$sysVar0)){
        sysVar0 <- sysVar1 <- NULL
        pAll <- interactions::cat_plot(model, pred=sysVar0, modx=sysVar1, y.label="Prob Profile = 2", x.label=sysVar0name, legend.main=sysVar1name, colors="Greys", interval=T)
      } else {
        sysVar0 <- sysVar1 <- NULL
        pAll <- interactions::interact_plot(model, pred=sysVar0, modx=sysVar1, y.label="Prob Profile = 2", x.label=sysVar0name, legend.main=sysVar1name, colors="Greys", interval=T)
      }
    }
    
    if(is.factor(basedata$sysVar0)){
      if(n_profiles == 3){
        ## 3 profiles categorical
        sysVar0 <- levels(basedata$sysVar0)
        sysVar1 <- levels(basedata$sysVar1)
        
        temp <- expand.grid(sysVar0, sysVar1)
        colnames(temp) <- c("sysVar0", "sysVar1")
        
        prob <- stats::predict(model, newdata=temp, "probs")
        colnames(prob) <- c("P1","P2","P3")
        temp2 <- cbind(prob, temp)
        
        vars1 <- c("P1", "sysVar0", "sysVar1")
        prob1 <- temp2[vars1]
        colnames(prob1) <- c("P1", sysVar0name, sysVar1name)
        
        vars2 <- c("P2", "sysVar0", "sysVar1")
        prob2 <- temp2[vars2]
        colnames(prob2) <- c("P2", sysVar0name, sysVar1name)
        
        vars3 <- c("P3", "sysVar0", "sysVar1")
        prob3 <- temp2[vars3]
        colnames(prob3) <- c("P3", sysVar0name, sysVar1name)
        
        pAll <- list()
        pAll[[1]] <- ggplot(data=prob1, aes_string(x=sysVar0name, y="P1", fill=sysVar1name)) +
          geom_bar(stat="identity", position=position_dodge()) +
          scale_fill_grey() +
          labs(title="Profile-1", y="Probabilty", x=sysVar0name) +
          scale_linetype_discrete(name=sysVar1name)
        
        pAll[[2]]<- ggplot(data=prob2, aes_string(x=sysVar0name, y="P2", fill=sysVar1name)) +
          geom_bar(stat="identity", position=position_dodge()) +
          scale_fill_grey() +
          labs(title="Profile-2", y="Probabilty", x=sysVar0name) +
          scale_linetype_discrete(name=sysVar1name)
        
        pAll[[3]] <- ggplot(data=prob3, aes_string(x=sysVar0name, y="P3", fill=sysVar1name)) +
          geom_bar(stat="identity", position=position_dodge()) +
          scale_fill_grey() +
          labs(title="Profile-3", y="Probabilty", x=sysVar0name) +
          scale_linetype_discrete(name=sysVar1name)
      }
      
      if(n_profiles == 4){
        ## 4 profiles categorical
        sysVar0 <- levels(basedata$sysVar0)
        sysVar1 <- levels(basedata$sysVar1)
        
        temp <- expand.grid(sysVar0, sysVar1)
        colnames(temp) <- c("sysVar0", "sysVar1")
        
        prob <- stats::predict(model, newdata=temp, "probs")
        colnames(prob) <- c("P1","P2","P3","P4")
        temp2 <- cbind(prob, temp)
        
        vars1 <- c("P1", "sysVar0", "sysVar1")
        prob1 <- temp2[vars1]
        colnames(prob1) <- c("P1", sysVar0name, sysVar1name)
        
        vars2 <- c("P2", "sysVar0", "sysVar1")
        prob2 <- temp2[vars2]
        colnames(prob2) <- c("P2", sysVar0name, sysVar1name)
        
        vars3 <- c("P3", "sysVar0", "sysVar1")
        prob3 <- temp2[vars3]
        colnames(prob3) <- c("P3", sysVar0name, sysVar1name)
        
        vars4 <- c("P4", "sysVar0", "sysVar1")
        prob4 <- temp2[vars4]
        colnames(prob4) <- c("P4", sysVar0name, sysVar1name)
        
        pAll <- list()
        pAll[[1]] <- ggplot(data=prob1, aes_string(x=sysVar0name, y="P1", fill=sysVar1name)) +
          geom_bar(stat="identity", position=position_dodge()) +
          scale_fill_grey() +
          labs(title="Profile-1", y="Probabilty", x=sysVar0name) +
          scale_linetype_discrete(name=sysVar1name)
        
        pAll[[2]] <- ggplot(data=prob2, aes_string(x=sysVar0name, y="P2", fill=sysVar1name)) +
          geom_bar(stat="identity", position=position_dodge()) +
          scale_fill_grey() +
          labs(title="Profile-2", y="Probabilty", x=sysVar0name) +
          scale_linetype_discrete(name=sysVar1name)
        
        pAll[[3]] <- ggplot(data=prob3, aes_string(x=sysVar0name, y="P3", fill=sysVar1name)) +
          geom_bar(stat="identity", position=position_dodge()) +
          scale_fill_grey() +
          labs(title="Profile-3", y="Probabilty", x=sysVar0name) +
          scale_linetype_discrete(name=sysVar1name)
        
        pAll[[4]] <- ggplot(data=prob3, aes_string(x=sysVar0name, y="P4", fill=sysVar1name)) +
          geom_bar(stat="identity", position=position_dodge()) +
          scale_fill_grey() +
          labs(title="Profile-4", y="Probabilty", x=sysVar0name) +
          scale_linetype_discrete(name=sysVar1name)
      }
    }
    
    if(is.numeric(basedata$sysVar0)){
      ### >2 continuous
      sysVar0L <- mean(basedata$sysVar0, na.rm=T) - stats::sd(basedata$sysVar0, na.rm=T)
      sysVar0H <- mean(basedata$sysVar0, na.rm=T) + stats::sd(basedata$sysVar0, na.rm=T)
      sysVar1L <- mean(basedata$sysVar1, na.rm=T) - stats::sd(basedata$sysVar1, na.rm=T)
      sysVar1H <- mean(basedata$sysVar1, na.rm=T) + stats::sd(basedata$sysVar1, na.rm=T)
      
      dataTemp<- matrix(c(sysVar0L, sysVar0H, sysVar0L, sysVar0H, sysVar1L, sysVar1L, sysVar1H, sysVar1H), nrow=4, ncol=2)
      dataTemp2 <- data.frame(dataTemp)
      colnames(dataTemp2) <- c("sysVar0", "sysVar1")
      
      prob <- stats::predict(model, newdata=dataTemp2, "probs")
      prob <- data.frame(stats::predict(model, newdata=dataTemp2, "probs"))
      prob$sysVar0 <- c(1,2,1,2)
      prob$sysVar1 <- c(1,1,2,2)
      prob$sysVar0 <- factor(prob$sysVar0, levels=c(1,2), labels=c("Low", "High"))
      prob$sysVar1 <- factor(prob$sysVar1, levels=c(1,2), labels=c("Low", "High"))
      
      if(n_profiles==3){
        vars1 <- c("X0", "sysVar0", "sysVar1")
        prob1 <- prob[vars1]
        vars2 <- c("X1", "sysVar0", "sysVar1")
        prob2 <- prob[vars2]
        vars3 <- c("X2", "sysVar0", "sysVar1")
        prob3 <- prob[vars3]
        
        pAll <- list()
        pAll[[1]] <- ggplot(prob1, aes_string(x = "sysVar0", y="X0", group="sysVar1")) + 
          geom_line(aes_string(linetype="sysVar1")) +
          ylim(0,1) +
          labs(title="Profile-1",y="Probabilty", x=sysVar0name) + 
          scale_linetype_discrete(name=sysVar1name)
        
        pAll[[2]] <- ggplot(prob2, aes_string(x = "sysVar0", y="X1", group="sysVar1")) + 
          geom_line(aes_string(linetype="sysVar1")) +
          ylim(0,1) +
          labs(title="Profile-2",y="Probabilty", x=sysVar0name) + 
          scale_linetype_discrete(name=sysVar1name)
        
        pAll[[3]] <- ggplot(prob3, aes_string(x = "sysVar0", y="X2", group="sysVar1")) + 
          geom_line(aes_string(linetype="sysVar1")) +
          ylim(0,1) +
          labs(title="Profile-3",y="Probabilty", x=sysVar0name) + 
          scale_linetype_discrete(name=sysVar1name)
      }
      
      if(n_profiles==4){
        vars1 <- c("X0", "sysVar0", "sysVar1")
        prob1 <- prob[vars1]
        vars2 <- c("X1", "sysVar0", "sysVar1")
        prob2 <- prob[vars2]
        vars3 <- c("X2", "sysVar0", "sysVar1")
        prob3 <- prob[vars3]
        vars4 <- c("X3", "sysVar0", "sysVar1")
        prob4 <- prob[vars4]
        
        pAll <- list()
        pAll[[1]] <- ggplot(prob1, aes_string(x = "sysVar0", y="X0", group="sysVar1")) + 
          geom_line(aes_string(linetype="sysVar1")) +
          ylim(0,1) +
          labs(title="Profile-1",y="Probabilty", x=sysVar0name) + 
          scale_linetype_discrete(name=sysVar1name)
        
        pAll[[2]]  <- ggplot(prob2, aes_string(x = "sysVar0", y="X1", group="sysVar1")) + 
          geom_line(aes_string(linetype="sysVar1")) +
          ylim(0,1) +
          labs(title="Profile-2",y="Probabilty", x=sysVar0name) + 
          scale_linetype_discrete(name=sysVar1name)
        
        pAll[[3]]  <- ggplot(prob3, aes_string(x = "sysVar0", y="X2", group="sysVar1")) + 
          geom_line(aes_string(linetype="sysVar1")) +
          ylim(0,1) +
          labs(title="Profile-3",y="Probabilty", x=sysVar0name) + 
          scale_linetype_discrete(name=sysVar1name)
        
        pAll[[4]]  <- ggplot(prob4, aes_string(x = "sysVar0", y="X3", group="sysVar1")) + 
          geom_line(aes_string(linetype="sysVar1")) +
          ylim(0,1) +
          labs(title="Profile-4",y="Probabilty", x=sysVar0name) + 
          scale_linetype_discrete(name=sysVar1name)
      }
    } 
  }
  
  if(printPlots==T){print(pAll)}  
  return(pAll)
}	




