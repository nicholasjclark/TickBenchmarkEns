
<!-- README.md is generated from README.Rmd. Please edit that file -->
# *TickBenchmarkEnsemble*

A Bayesian log-Gamma weighted ensemble forecast of *Ambloyomma americanum* larval density for the [NEON EFI Forecasting Challenge 2022](https://projects.ecoforecast.org/neon4cast-docs/).

## Models

The method is described as follows: First a set of benchmark models are fitted to each site-level time series; these models consist of a robust STL decomposition with AR errors, a TBATS model and a GAM (Gamma response distribution with smooth additive predictors of both the location and the scale). Second, the fitted and point forecast predicted values from each model are used as predictors in a weighted ensemble Bayesian log-Gamma model in which the weights are estimated as multivariate normal regression coefficients alongside a random walk process to capture any remaining autocorrelation. Parameters of this model are estimated using MCMC sampling via a Gibbs sampler in the software `JAGS`.

The entire modelling process is outlined in the `R/benchmark_ens_jags.R` script

## Purpose

There are limitations with these models in that (1) missing values must first be imputed prior to fitting the STL and TBATS models, and this imputation scheme will have major impacts on the resulting forecasts; (2) the current modelling framework does not use partial pooling among locations and; (3) no environmental covariates are included. The purpose of the modelling framework is to identify useful benchmark models that can perform reasonably well in forecasting tick densities so that more detailed population-based or mechanistic models that incorporate biological realism can be tested and improved.

## License

None whatsoever. Please use, adapt and distribute this code widely so that the models can be improved and our collective knowledge of tick abundance distributions can accelerate.
