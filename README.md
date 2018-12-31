# rties
## Tools for Modeling Interpersonal Dynamics
This package grew out of my research on temporal interpersonal emotion systems (TIES; See TIES_PSPB_2011.pdf in the documentation folder), hence the name "rties". It provides tools for using a suite of models to investigate temporal processes in bivariate (e.g., dyadic) systems. Although my focus is on interpersonal emotional processes, the models could be used for any bivariate system. The general approach is to model - one dyad at a time - the dynamics of a variable that is assessed repeatedly from both partners, extract the parameter estimates for each dyad, and then use those parameter estimates to either predict, or be predicted by, another variable of interest. Currently, 3 models are supported: 1) inertia-coordination, 2) patterned-slopes, and 3) a coupled-oscillator. The models are currently restricted to distinguishable dyads and, if using the dynamics to predict outcomes, those outcome variables must be normally distributed.
Documentation is available in the form of pdf documents in the documentation folder.

## Vignettes available at RPubs
overview_data_prep: https://rpubs.com/ebmtnprof/455051
