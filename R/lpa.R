

#################### inspectProfiles

#' A wrapper for tidyLPA that provides information to help decide how many profiles to use for  subsequent rties analyses. 
#'
#' The function prints out the number of dyads in each profile for a specified number of profiles. It also prints out the tidyLPA print_profiles figure and prototypical model-predicted trajectories for each profile. For the inertia-coordination model, it produces sets of prototypical examples by using the inertia-coordination parameters to predict temporal trajectories, with random noise added at each temporal step. This process is required because the inertia-coordination model only represents local dynamics and predictions bear no resemblance to observed variables without the addition of noise. For the coupled-oscillator, this step is not necessary and one prototypical trajectory is plotted for each profile.
#' 
#' @param whichModel The name of the model that is being investigated (e.g., "inertCoord" or "clo")
#' @param prepData A dataframe that was produced with the "dataPrep" function.
#' @param paramEst A dataframe created by either indivInertCoord or indivClo containing the parameter estimates for each dyad.
#' @param hasNA Indicates whether paramEst contains missing data (TRUE) or not (FALSE). If it does (hasNA=T) then random-forest imputation is done prior to estimating the latent profiles.
#' @param n_profiles The number of latent profiles.
#' @param dist0name An optional name for the level-0 of the distinguishing variable (e.g., "Women"). Default is dist0.
#' @param dist1name An optional name for the level-1 of the distinguishing variable (e.g., "Men"). Default is dist1
#' @param obsName An optional name for the state variables being plotted (e.g., "heart rate"). Default is obsName.
#' @param minMax An optional vector with desired minimum and maximum quantiles to be used for setting the y-axis range on the plots, e.g., minMax <- c(.1, .9) would set the y-axis limits to the 10th and 90th percentiles of the observed state variables. If not provided, the default is to use the minimum and maximum observed values of the state variables.
#' 
#' @return The function prints: 1) the number of dyads in each profile, 2) the figure created by tidyLPA's "plot_profiles", which gives a visual indicator of the quality of the LPA solution, and 3) prototypical model-predicted trajectories for each profile. 

#' @import ggplot2
#' @export

inspectProfiles <- function(whichModel, prepData, paramEst, hasNA, n_profiles, dist0name=NULL, dist1name=NULL, obsName=NULL, minMax=NULL)
{   
   if(whichModel == "clo" & hasNA == TRUE){
    lpa <- paramEst %>%
        select(obs_0:p_d1_1) %>%
        tidyLPA::single_imputation(method="missForest") %>%
        tidyLPA::estimate_profiles(n_profiles)
  } else if (whichModel == "clo" & hasNA == FALSE){
  	  lpa <- paramEst %>%
        select(obs_0:p_d1_1) %>%
        tidyLPA::estimate_profiles(n_profiles)
    } else if (whichModel == "inertCoord" & hasNA == TRUE){
  	    lpa <- paramEst %>%
          select(inert1, coord1, coord0, inert0) %>%
          tidyLPA::single_imputation(method="missForest") %>%
          tidyLPA::estimate_profiles(n_profiles)
       } else if (whichModel == "inertCoord" & hasNA == FALSE){
  	      lpa <- paramEst %>%
            select(inert1, coord1, coord0, inert0) %>%
            tidyLPA::estimate_profiles(n_profiles)
         } else 
        print("Model must be inertCoord or clo and hasNA must be T or F")      
        
  profileData <- as.data.frame(tidyLPA::get_data(lpa))
  temp1 <- profileData[!duplicated(profileData$id), ]
  nPerProfile <- data.frame(table(temp1$Class))
  colnames(nPerProfile) <- c("Profile", "Frequency")
  print(nPerProfile)  
  
  tidyLPA::plot_profiles(lpa)
   
  if(whichModel=="clo") {
  cloPlotTrajInternal(prepData=prepData, paramEst=paramEst, n_profiles=n_profiles,dist0name=dist0name, dist1name=dist1name, obsName=obsName, minMax=minMax)
  } else if (whichModel=="inertCoord"){
  	  inertCoordPlotTrajInternal(prepData=prepData, paramEst=paramEst, n_profiles=n_profiles, dist0name=dist0name, dist1name=dist1name, obsName=obsName, minMax=minMax)
  } else
  print("Model must be inertCoord or clo")
  
  return(profileData)
}

