
######## This file includes the function needed to estimate the latent profiles based on IC or CLO model parameters


#################### inspectProfiles

#' Provides information to help decide how many profiles to use for subsequent rties analyses. 
#'
#' The function prints out the number of dyads in each profile for a specified number of profiles. It also prints out: 1) a figure showing the best clustering solution as indicated by BIC (e.g., the observed data separated into clusters, produced by mclust), 2) a line plot showing the content of the best solution (e.g., the mean parameter estimates for each profile) and 3) prototypical model-predicted trajectories for each profile. For the inertia-coordination model, it produces sets of prototypical examples by using the inertia-coordination parameters to predict temporal trajectories, with random noise added at each temporal step. This process is required because the inertia-coordination model only represents local dynamics and predictions bear no resemblance to observed variables without the addition of noise. An optional argument, "seed" sets the seed for the random number generator, so you can get the same plots each time. If the "seed" argument is used, then only one plot per profile is produced. For the coupled-oscillator, this step is not necessary and one prototypical trajectory is plotted for each profile.
#' 
#' @param whichModel The name of the model that is being investigated (e.g., "inertCoord" or "clo")
#' @param time_lag The time-lag used for the inertia-coordination model. This is required for the "inertCoord" model and not relevant to the "clo" model.
#' @param prepData A dataframe that was produced with the "dataPrep" function.
#' @param paramEst A dataframe created by either indivInertCoord or indivClo containing the parameter estimates for each dyad.
#' @param n_profiles The number of latent profiles.
#' @param dist0name An optional name for the level-0 of the distinguishing variable (e.g., "Women"). Default is dist0.
#' @param dist1name An optional name for the level-1 of the distinguishing variable (e.g., "Men"). Default is dist1
#' @param minMax An optional vector with desired minimum and maximum quantiles to be used for setting the y-axis range on the plots, e.g., minMax <- c(.1, .9) would set the y-axis limits to the 10th and 90th percentiles of the observed state variables. If not provided, the default is to use the minimum and maximum observed values of the state variables.
#' @param set.seed An optional integer argument that sets the seed of R's random number generator to create reproducible trajectories. If used, the number of plots produced is set to one for each profile.
#' 
#' @return The function returns a dataframe called "profileData" that contains the profile classification for each dyad. 

#' @import ggplot2
#' @import mclust
#' @export

inspectProfiles <- function(whichModel, prepData, paramEst, n_profiles, dist0name=NULL, dist1name=NULL, minMax=NULL, seed = NULL)
{  
	profileData <- list() 
	paramEst <- paramEst[complete.cases(paramEst), ]
   
   if(whichModel == "clo"){
  	  params <- subset(paramEst, select=c(obs_0:p_d1_1))
  	  lpa <- Mclust(params, G=n_profiles)
    } else if (whichModel == "inertCoord"){
  	    params <- subset(paramEst, select=c(inert1, coord1, coord0, inert0))
  	    lpa <- Mclust(params, G=n_profiles)
         } else 
        print("Model must be inertCoord or clo") 
	  
	  if(! is.null(seed)) {set.seed = seed}
   
   # profileData
   profileData$profile <- factor(lpa$classification)
   profileData$profileN <- as.numeric(lpa$classification) - 1
   profileData$dyad <- paramEst$dyad
   profileData <- as.data.frame(profileData)
   
   # classification table     
  print(table(lpa$classification)) 
  
  # plot quality of solution
  dr <- MclustDR(lpa, lambda=1)
  plot(dr, what ="contour")
  
  # plot content of solution
  means <- as.data.frame(lpa$parameters$mean)
  means$varNames <- rownames(means)
  means$var <- c(1:dim(means)[1])
  meansL <- reshape(means, idvar="varNames", varying=list(1:n_profiles), timevar="profile", sep="", direction="long")

  print(ggplot(data=meansL, aes(x=varNames, y=V1, group=profile)) +
		geom_line(aes(colour=as.factor(profile))))

  if(whichModel=="clo") {
  cloPlotTrajInternal(prepData=prepData, paramEst=paramEst, n_profiles=n_profiles, dist0name=dist0name, dist1name=dist1name, minMax=minMax, seed=seed)
  } else if (whichModel=="inertCoord"){
  	  inertCoordPlotTrajInternal(prepData=prepData, paramEst=paramEst, n_profiles=n_profiles, dist0name=dist0name, dist1name=dist1name, minMax=minMax, seed=seed)
  } else
  print("Model must be inertCoord or clo")
  
  return(profileData)
}

