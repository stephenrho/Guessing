# Materials and Analysis Code for Rhodes, Cowan, Hardman & Logie: Informed Guessing in Change Detection

This repository contains the materials (`materials` folder), data, and analysis code (`analysis` folder) for the manuscript *Informed Guessing in Change Detecton*.

The scripts in the `analysis` folder do the following:

* `createDataLists.R` restructures the data into lists for use with R2jags.
* `models.R` creates the JAGS models and writes to a folder `models`.
* `samlplesExp[1-4].R` uses [JAGS](http://mcmc-jags.sourceforge.net/) to sample from the joint posterior of each model under consideration.
* `PP_check.R` calculates the overlap between the posterior predictive distribution of the informed model and the raw data.
* `guessingFunctions.R` contains miscellaneous functions used in restructuring posterior distributions/ plotting. Some functions are from ([Kruschke, 2015](http://www.indiana.edu/~kruschke/DoingBayesianDataAnalysis/)).
* `conventional_k.R` conducts a standard analysis of *k* using the models of [Pashler (1988)](https://www.ncbi.nlm.nih.gov/pubmed/3226885) and [Cowan et al. (2013)](https://www.ncbi.nlm.nih.gov/pubmed/22905929).
* `relativeLuminance.R` calculates the difference in average luminance between study and test arrays and performs logistic mixed effects analysis to assess role of luminance in driving the probability of responding change on a given trial for each experiment (response to reviewer comment).
* `variableU.R` examines the posterior of a model that estimated a separate guessing rate for each set size (response to reviewer). Plots of the resulting estimates for each experiment can be found in the folder `vary_u_figs` with posterior samples from the guessing rate from the informed guessing model for comparison.
* `individualK.R` can be used to look at participant level estimates of *k*.
* `waicFunctions.R` and `waic.R` contain code written to calculate [WAIC](http://www.jmlr.org/papers/v11/watanabe10a.html) with the posterior chains of each model.
