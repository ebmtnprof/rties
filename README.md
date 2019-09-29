# rties

## Please cite as : 
Butler, E. A. & Barnard, J. K. (in press). Quantifying interpersonal dynamics for studying socio-emotional processes and adverse health behaviors. Psychosomatic Medicine.
-- see documentation folder for the pdf

## Tools for Modeling Interpersonal Dynamics
This package grew out of our research on temporal interpersonal emotion systems (TIES; See TIES_PSPB_2011.pdf in the documentation folder), hence the name "rties". It provides tools for using a suite of models to explore temporal processes in bivariate (e.g., dyadic) systems. Although our focus is on interpersonal emotional processes, the models could be used for any bivariate system (e.g., heart rate and skin conductance measured over time within individuals). The general approach is to model - one dyad at a time - the dynamics of a variable that is assessed repeatedly from both partners, extract the parameter estimates for each dyad, and then use those parameter estimates to either predict, or be predicted by, another variable of interest. Currently, two models are supported: 1) inertia-coordination, and 2) a coupled-oscillator. Documentation is available in the form of pdfs in the "documentation" folder (relevant published manuscripts) and vignettes published as RPubs (detailed "how to" vignettes listed below).

Note: We are updating the package regularly, fixing bugs and adding functionality, so please check back often for updates and email us if you have problems.

## Vignettes available at RPubs

These vignettes provide an overview of rties functionality and work flow. Many of the functions have optional arguments that are not documentd in the vignettes. For full information, use ?functionName (where functionName is the name of the function you want information for).

overview_data_prep_V03: http://rpubs.com/ebmtnprof/484725

inertia_coordination_V03: http://rpubs.com/ebmtnprof/495166

coupled_oscillator_V03: http://rpubs.com/ebmtnprof/495181
