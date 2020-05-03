
######## This file includes the function needed to estimate the latent profiles based on IC or CLO model parameters


#################### inspectProfiles

#' Provides information to help decide how many profiles to use for subsequent rties analyses. 
#'
#' The function prints out the number of dyads in each profile for a specified number of profiles. It also prints out: 1) a figure showing the best clustering solution as indicated by BIC (e.g., the observed data separated into clusters, produced by mclust), 2) a line plot showing the content of the best solution (e.g., the mean parameter estimates for each profile) and 3) prototypical model-predicted trajectories for each profile. For the inertia-coordination model, it produces sets of prototypical examples by using the inertia-coordination parameters to predict temporal trajectories, with random noise added at each temporal step. This process is required because the inertia-coordination model only represents local dynamics and predictions bear no resemblance to observed variables without the addition of noise. An optional argument, "seed" sets the seed for the random number generator, so you can get the same plots each time. If the "seed" argument is used, then only one plot per profile is produced. For the coupled-oscillator, this step is not necessary and one prototypical trajectory is plotted for each profile.
#' 
#' @param whichModel The name of the model that is being investigated (e.g., "inertCoord" or "clo")
#' @param prepData A dataframe that was produced with the "dataPrep" function.
#' @param paramEst A dataframe created by either indivInertCoord or indivClo containing the parameter estimates for each dyad.
#' @param n_profiles The number of latent profiles.
#' @param dist0name An optional name for the level-0 of the distinguishing variable (e.g., "Women"). Default is dist0.
#' @param dist1name An optional name for the level-1 of the distinguishing variable (e.g., "Men"). Default is dist1
#' @param plot_obs_name An optional name for the observed state variable to appear on plots (e.g., "Emotional Experience").
#' @param minMax An optional vector with desired minimum and maximum quantiles to be used for setting the y-axis range on the plots, e.g., minMax <- c(.1, .9) would set the y-axis limits to the 10th and 90th percentiles of the observed state variables. If not provided, the default is to use the minimum and maximum observed values of the state variables.
#' @param time_length An optional value specifying how many time points to plot across. Default is the 75th percentile for the time variable.
#' @param numPlots Only relevant for the inertCoord model. An optional value controlling how many random examples of each profile are produced. Default is 3.
#' @param seed Only relevant for the inertCoord model. An optional integer argument that sets the seed of R's random number generator to create reproducible trajectories. If used, the "numPlots" can be set to one - otherwise each plot is replicated 3 times.
#' @examples
#' data <- rties_ExampleDataShort
#' newData <- dataPrep(basedata=data, dyadId="couple", personId="person", 
#'  obs_name="dial", dist_name="female", time_name="time")
#' taus <-c(2,3)
#' embeds <- c(3,4)
#' delta <- 1
#' derivs <- estDerivs(prepData=newData, taus=taus, embeds=embeds, delta=delta, idConvention=500)
#' clo <- indivClo(derivData=derivs$data, whichModel="coupled")
#' profiles <- inspectProfiles(whichModel="clo", prepData=newData, paramEst=clo$params, n_profiles=2)
#' head(profiles)
#' 
#' @return A dataframe called "profileData" that contains the profile classification for each dyad. 

#' @import ggplot2
#' @import mclust

#' @export

inspectProfiles <- function(whichModel, prepData, paramEst, n_profiles, dist0name=NULL, dist1name=NULL, plot_obs_name = NULL, minMax=NULL, time_length=NULL, numPlots=NULL, seed = NULL)
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
  if(is.null(numPlots)) {numPlots <- 3}
  if(!is.null(seed)) {seed = seed}
  
  profileData <- list() 
	paramEst <- paramEst[stats::complete.cases(paramEst), ]
   
   if(whichModel == "clo"){
  	  vars1 <- c("obs_0","d1_0","p_obs_0","p_d1_0","obs_1","d1_1","p_obs_1","p_d1_1")
  	  params <- paramEst[vars1]
  	  lpa <- Mclust(params, G=n_profiles)
    } else if (whichModel == "inertCoord"){
  	    vars2 <- c("inert1", "coord1", "coord0", "inert0")
  	    params <- paramEst[vars2]
  	    lpa <- Mclust(params, G=n_profiles)
         } else 
        print("Model must be inertCoord or clo") 
   
   # profileData
   profileData$profile <- factor(lpa$classification)
   profileData$profileN <- as.numeric(lpa$classification) - 1
   profileData$dyad <- paramEst$dyad
   profileData <- as.data.frame(profileData)
   
   # classification table     
  print(table(lpa$classification)) 
  
  # plot quality of solution
  dr <- mclust::MclustDR(lpa, lambda=1)
  graphics::plot(dr, what ="contour")
  
  # plot content of solution
  means <- as.data.frame(lpa$parameters$mean)
  means$varNames <- rownames(means)
  means$var <- c(1:dim(means)[1])
  meansL <- stats::reshape(means, idvar="varNames", varying=list(1:n_profiles), timevar="profile", sep="", direction="long")

  profile <- NULL
  print(ggplot(data=meansL, aes_string(x="varNames", y="V1", group="profile")) +
		geom_line(aes(colour=as.factor(profile))))

  if(whichModel=="clo") {
  cloPlotTraj(prepData=prepData, paramEst=paramEst, n_profiles=n_profiles, dist0name=dist0name, dist1name=dist1name, plot_obs_name=plot_obs_name, minMax=minMax, time_length=time_length)
  } else if (whichModel=="inertCoord"){
  	  inertCoordPlotTraj(prepData=prepData, paramEst=paramEst, n_profiles=n_profiles, dist0name=dist0name, dist1name=dist1name, plot_obs_name=plot_obs_name, minMax=minMax, time_length=time_length, numPlots=numPlots, seed=seed)
  } else
  print("Model must be inertCoord or clo")
  
  return(profileData)
}
