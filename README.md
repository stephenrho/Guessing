# Materials and Analysis Code for Rhodes, Cowan, Hardman & Logie: Informed Guessing in Change Detection

This repository contains the materials (`materials` folder), data, and analysis code (`analysis` folder) for the manuscript *Informed Guessing in Change Detecton*.

The scripts in the `analysis` folder do the following:

* `models.R` creates the JAGS models and writes to a folder `models`
* `samlplesExp[1-4].R` uses [JAGS](http://mcmc-jags.sourceforge.net/) to sample from the joint posterior of each model under consideration.
* `PP_check.R` calculates the overlap between the posterior predictive distribution and the raw data.
* `guessingFunctions.R` contains miscellaneous functions used in restructuring posterior distributions/ plotting. Some functions are from ([Kruschke, 2015](http://www.indiana.edu/~kruschke/DoingBayesianDataAnalysis/)).
* `conventional_k.R` conducts a standard analysis of *k* using the models of [Pashler (1988)](https://www.ncbi.nlm.nih.gov/pubmed/3226885) and [Cowan et al. (2013)](https://www.ncbi.nlm.nih.gov/pubmed/22905929).
