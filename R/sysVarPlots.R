
############### sysVarInPlots

#' Produces plots for interpreting the results from sysVarIn.
#'
#' @param fullData A dataframe created by the "makeFullData" function.
#' @param sysVar_name The name of the variable in the dataframe that contains the system variable.
#' @param sysVarType Whether the system variable is "dyadic", which means both partners have the same score, or "indiv" which means the partners can have different scores
#' @param n_profiles The number of latent profiles.
#' @param testModel The name of the model that is being interpreted (e.g., sysIn$models$sysVarInteract). Only needed when the system variable is "indiv" (e.g., individual scores for each partner)
#' @param dist0name An optional name for the level-0 of the distinguishing variable (e.g., "Women"). Default is dist0.
#' @param dist1name An optional name for the level-1 of the distinguishing variable (e.g., "Men"). Default is dist1
#' @param printPlots If true (the default) plots are displayed on the screen.
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
#' sysVarInPlots(fullData=fullData, sysVar_name="conflict", sysVarType="indiv", 
#' n_profiles=2, testModel=sysIn$models$sysVarInteract)
#' 
#' @return Single plots or a list of plots (depending on the model that is being interpreted).

#' @import ggplot2
#' @export

sysVarInPlots <- function(fullData, sysVar_name, sysVarType, n_profiles, testModel=NULL, dist0name=NULL, dist1name=NULL, printPlots=T){
  basedata <- fullData
  basedata <- basedata[stats::complete.cases(basedata), ]
  colnames(basedata)[colnames(basedata)== sysVar_name] <- "sysVar" 	
  if(n_profiles > 4) {message("plots are not provided if there are more than 4 profiles") }
  
  if(is.null(dist0name)){dist0name <- "dist0"}
  if(is.null(dist1name)){dist1name <- "dist1"}

  if(sysVarType == "dyadic"){
    pAll <- dyadic(basedata, sysVar_name)
  }

  if(sysVarType == "indiv"){
    vars1 <- c("dyad", "sysVar", "dist0", "profileN")
    data1 <- basedata[vars1]
    data2 <-  stats::reshape(data1, idvar="dyad", timevar = "dist0", direction= "wide")   
    dyad <- sysVar.1 <- profileN.1 <- sysVar.0 <- profileN.0 <- NULL
    data3 <- dplyr::rename(data2, dyad=dyad, sysVar0=sysVar.1, profileN1= profileN.1, sysVar1= sysVar.0, profileN=profileN.0)
    basedata <- data3[stats::complete.cases(data3), ]
    sysVar0name <- paste(sysVar_name, dist0name, sep="_")
    sysVar1name <- paste(sysVar_name, dist1name, sep="_")
    
    if(is.factor(basedata$sysVar0)){
      if(n_profiles == 2){
        pAll <- indiv2profilesCat(testModel, sysVar0name, sysVar1name)
       } 
      if(n_profiles == 3){
        pAll <- indiv3profilesCat(basedata, testModel, sysVar0name, sysVar1name)
        }
      if(n_profiles == 4){
        pAll <- indiv4profilesCat(basedata, testModel, sysVar0name, sysVar1name)
      }
    }
      
    if(is.numeric(basedata$sysVar0)){
      
      if(n_profiles == 2){
        pAll <- indiv2profilesCont(testModel, sysVar0name, sysVar1name)  
      }
      
      if(n_profiles > 2){
      sysVar0L <- mean(basedata$sysVar0, na.rm=T) - stats::sd(basedata$sysVar0, na.rm=T)
      sysVar0H <- mean(basedata$sysVar0, na.rm=T) + stats::sd(basedata$sysVar0, na.rm=T)
      sysVar1L <- mean(basedata$sysVar1, na.rm=T) - stats::sd(basedata$sysVar1, na.rm=T)
      sysVar1H <- mean(basedata$sysVar1, na.rm=T) + stats::sd(basedata$sysVar1, na.rm=T)
      dataTemp<- matrix(c(sysVar0L, sysVar0H, sysVar0L, sysVar0H, sysVar1L, sysVar1L, sysVar1H, sysVar1H), nrow=4, ncol=2)
      dataTemp2 <- data.frame(dataTemp)
      colnames(dataTemp2) <- c("sysVar0", "sysVar1")
      prob <- data.frame(stats::predict(testModel, newdata=dataTemp2, type="probs"))
      prob$sysVar0 <- c(1,2,1,2)
      prob$sysVar1 <- c(1,1,2,2)
      prob$sysVar0 <- factor(prob$sysVar0, levels=c(1,2), labels=c("Low", "High"))
      prob$sysVar1 <- factor(prob$sysVar1, levels=c(1,2), labels=c("Low", "High"))
      }
      
      if(n_profiles == 3){
        pAll <- indiv3profilesCont(prob, sysVar0name, sysVar1name)
        }
      if(n_profiles == 4){
        pAll <- indiv4profilesCont(prob, sysVar0name, sysVar1name)
        }
      }
    }
  if(printPlots==T){print(pAll)}  
  return(pAll)
}

########################## The following functions are called by sysVarInPlots

####### dyadic

#' Produces plots for sysVarIn when sysVar is dyadic.
#'
#' @param basedata A dataframe created internally by the "sysVarInPlots" function.
#' @param sysVar_name The name of the variable in the dataframe that contains the system variable.
#' 
#' @return A plot with the profiles on the y-axis and the system variable on the x-axis

dyadic <- function(basedata, sysVar_name){
  basedata <- basedata[!duplicated(basedata$dyad), ]	
  
  sysVar <- ..prop.. <- profile <- profileN <- NULL
  
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
  return(pAll)
}

####### indiv2profilesCat

#' Produces plots for sysVarIn when sysVar is categorical and there are 2 profiles
#'
#' @param testModel The model object created by sysVarIn for the interaction model (e.g., sysVarInteract)
#' @param sysVar0name The name created by sysVarInPlots referring to the system variable for partner-0.
#' @param sysVar1name The name created by sysVarInPlots referring to the system variable for partner-1.
#' 
#' @return A plot produced by the interactions package.

indiv2profilesCat <- function(testModel, sysVar0name, sysVar1name){
  sysVar0 <- sysVar1 <- NULL
  pAll <- interactions::cat_plot(testModel, pred=sysVar0, modx=sysVar1, y.label="Prob Profile = 2", x.label=sysVar0name, legend.main=sysVar1name, colors="Greys", interval=T)
  return(pAll)
} 

####### indiv2profilesCont

#' Produces plots for sysVarIn when sysVar is continuous and there are 2 profiles
#'
#' @param testModel The model object created by sysVarIn for the interaction model (e.g., sysVarInteract)
#' @param sysVar0name The name created by sysVarInPlots referring to the system variable for partner-0.
#' @param sysVar1name The name created by sysVarInPlots referring to the system variable for partner-1.
#' 
#' @return A plot produced by the interactions package.

indiv2profilesCont <- function(testModel,sysVar0name, sysVar1name) {
  sysVar0 <- sysVar1 <- NULL
  pAll <- interactions::interact_plot(testModel, pred=sysVar0, modx=sysVar1, y.label="Prob Profile = 2", x.label=sysVar0name, legend.main=sysVar1name, colors="Greys", interval=T)
  return(pAll)
}

####### indiv3profilesCat

#' Produces plots for sysVarIn when sysVar is categorical and there are 3 profiles

#' @param basedata A dataframe created internally by the "sysVarInPlots" function.
#' @param testModel The model object created by sysVarIn for the interaction model (e.g., sysVarInteract)
#' @param sysVar0name The name created by sysVarInPlots referring to the system variable for partner-0.
#' @param sysVar1name The name created by sysVarInPlots referring to the system variable for partner-1.
#' 
#' @return A list of 3 plots showing the simple slopes for each of the profiles.

indiv3profilesCat <- function(basedata, testModel, sysVar0name, sysVar1name){
  sysVar0 <- levels(basedata$sysVar0)
  sysVar1 <- levels(basedata$sysVar1)
  
  temp <- expand.grid(sysVar0, sysVar1)
  colnames(temp) <- c("sysVar0", "sysVar1")
  
  prob <- stats::predict(testModel, newdata=temp, "probs")
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
  
  return(pAll)
}

####### indiv4profilesCat

#' Produces plots for sysVarIn when sysVar is categorical and there are 4 profiles

#' @param basedata A dataframe created internally by the "sysVarInPlots" function.
#' @param testModel The model object created by sysVarIn for the interaction model (e.g., sysVarInteract)
#' @param sysVar0name The name created by sysVarInPlots referring to the system variable for partner-0.
#' @param sysVar1name The name created by sysVarInPlots referring to the system variable for partner-1.
#' 
#' @return A list of 4 plots showing the simple slopes for each of the profiles.

indiv4profilesCat <- function(basedata, testModel, sysVar0name, sysVar1name){
  sysVar0 <- levels(basedata$sysVar0)
  sysVar1 <- levels(basedata$sysVar1)
  
  temp <- expand.grid(sysVar0, sysVar1)
  colnames(temp) <- c("sysVar0", "sysVar1")
  
  prob <- stats::predict(testModel, newdata=temp, "probs")
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
  
  pAll[[4]] <- ggplot(data=prob4, aes_string(x=sysVar0name, y="P4", fill=sysVar1name)) +
    geom_bar(stat="identity", position=position_dodge()) +
    scale_fill_grey() +
    labs(title="Profile-4", y="Probabilty", x=sysVar0name) +
    scale_linetype_discrete(name=sysVar1name)
  
  return(pAll)
}

####### indiv3profilesCont

#' Produces plots for sysVarIn when sysVar is continuous and there are 3 profiles

#' @param prob A dataframe created internally by the "sysVarInPlots" function.
#' @param sysVar0name The name created by sysVarInPlots referring to the system variable for partner-0.
#' @param sysVar1name The name created by sysVarInPlots referring to the system variable for partner-1.
#' 
#' @return A list of 3 plots showing the simple slopes for each of the profiles.
#' 
indiv3profilesCont <- function(prob, sysVar0name, sysVar1name){
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
  
  return(pAll)
} 

####### indiv4profilesCont

#' Produces plots for sysVarIn when sysVar is continuous and there are 4 profiles

#' @param prob A dataframe created internally by the "sysVarInPlots" function.
#' @param sysVar0name The name created by sysVarInPlots referring to the system variable for partner-0.
#' @param sysVar1name The name created by sysVarInPlots referring to the system variable for partner-1.
#' 
#' @return A list of 4 plots showing the simple slopes for each of the profiles.
#' 
indiv4profilesCont <- function(prob, sysVar0name, sysVar1name){
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
  
  pAll[[4]] <- ggplot(prob4, aes_string(x = "sysVar0", y="X3", group="sysVar1")) + 
    geom_line(aes_string(linetype="sysVar1")) +
    ylim(0,1) +
    labs(title="Profile-4",y="Probabilty", x=sysVar0name) + 
    scale_linetype_discrete(name=sysVar1name)
  
  return(pAll)
} 

#######################################

############### sysVarOutPlots

#' Produces plots for interpreting the results from sysVarIn.
#'
#' @param fullData A dataframe created by the "makeFullData" function.
#' @param sysVar_name The name of the variable in the dataframe that contains the system variable.
#' @param sysVarType Whether the system variable is "dyadic", which means both partners have the same score, or "indiv" which means the partners can have different scores
#' @param testModel The name of the model that is being interpreted (e.g., sysIn$models$sysVarInteract). 
#' @param dist0name An optional name for the level-0 of the distinguishing variable (e.g., "Women"). Default is dist0.
#' @param dist1name An optional name for the level-1 of the distinguishing variable (e.g., "Men"). Default is dist1
#' @param binomial Whether the system variable is binomial. Default is false.
#' @examples
#' # See vignettes for examples.
#' 
#' @return Single plots or a list of plots (depending on the model that is being interpreted).

#' @import ggplot2
#' @export

sysVarOutPlots <- function(fullData, sysVar_name, sysVarType, testModel, dist0name=NULL, dist1name=NULL, binomial=F){
  
  basedata <- fullData
  basedata <- basedata[stats::complete.cases(basedata), ]
  if(is.null(dist0name)){dist0name <- "dist0"}
  if(is.null(dist1name)){dist1name <- "dist1"}
  colnames(basedata)[colnames(basedata)== sysVar_name] <- "sysVar" 	
  basedata$dist <- factor(basedata$dist0, labels=c(dist1name, dist0name))
  
  profile <- sysVar <- ..prop.. <- NULL
  
  if(binomial==F){
    if(sysVarType == "dyadic")
    {
      resid <- data.frame(resid(testModel))
      colnames(resid) <- "resid"
      pAll <- list()
      
      pAll[[1]] <- ggplot(resid, aes(x=resid)) +
        geom_histogram(color="black", fill="grey")
      
      pAll[[2]] <- ggplot(basedata, aes(x=profile, y=sysVar)) +
        geom_boxplot() + 
        ylab(sysVar_name)
    }
    
    if(sysVarType == "indiv")
    {
      resid <- data.frame(resid(testModel))
      colnames(resid) <- "resid"
      pAll <- list()
      
      pAll[[1]] <- ggplot(resid, aes(x=resid)) +
        geom_histogram(color="black", fill="grey")
      
      temp <- sjPlot::plot_model(testModel, type="pred", terms=c("profile", "dist"), colors="gs", y.label=sysVar_name)
      pAll[[2]] <- temp + ylab(sysVar_name) 
    }
  }
  
  if(binomial==T){
    if(sysVarType == "dyadic")
    {
      label <- paste("Proportions", sysVar_name, "= 0 or 1 in each profile", sep=" ")
      pAll <- ggplot(basedata) +
        geom_bar(aes(x = profile, y = ..prop.., group = sysVar)) +
        facet_wrap(~ sysVar) +
        ylab(label)
    }
    
    if(sysVarType == "indiv")
    {
      temp <- sjPlot::plot_model(testModel, type="pred", terms=c("profile", "dist"), colors="gs", y.label=sysVar_name)
      pAll <- temp + ylab(sysVar_name) 
    }
  }
  return(pAll)
}




  