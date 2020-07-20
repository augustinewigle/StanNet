StanNet is a variant of the BUGSnet package (Bayesian inference Using Gibbs Sampling to conduct NETwork meta-analysis) which uses Stan in place of JAGS for MCMC sampling. StanNet allows the user to perform Bayesian Network Meta-Analysis and easily produce summaries of their data and results.
The package supports additional data exploration, simulation, and diagnostics not originally present in BUGSnet. However, it currently only supports two-arm trials and does not support meta-regression or inconsistency tests.
Users of this package must first download Stan in R on their computer prior to use. Instructions to do so can be found [here](https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started).