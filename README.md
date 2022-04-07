
<!-- README.md is generated from README.Rmd. Please edit that file -->
# *TickBenchmarkEnsemble*

A Bayesian log-Gamma weighted ensemble forecast of *Ambloyomma americanum* larval density for the [NEON EFI Forecasting Challenge 2022](https://projects.ecoforecast.org/neon4cast-docs/).

## Models

The method is described as follows: First a a Generalised Additive Model (GAM) with an assumed Tweedie response distribution is fitted to a training sample of each individual time series. This model estimates smooth functions for seasonality and long-term trend. Next, posterior expectation distributions for the training series are extracted from the GAM as a way of 'imputing' missing values with some respect of uncertainty around those values. Imputed in-sample training series on the log scale are then used to fit two additional benchmark sub-models to each site-level time series; these sub-models consist of a robust STL decomposition with AR errors and a TBATS model. Second, the fitted values and testing horizon point forecasts (on the log scale) from each model (i.e. the GAM, the STL-AR and the TBATS) are used as predictors in a weighted ensemble Bayesian log-Gamma model in which the weights are estimated as multivariate normal regression coefficients. Parameters of this model are estimated using MCMC sampling via a Gibbs sampler in the software `JAGS`. The estimated weight distributions are then used to update the forecast by training each of the above three models on the entire available time series for each site.

The entire modelling process is outlined in the `R/benchmark_ens_jags.R` script

## Purpose

There are limitations with these models in that (1) missing values must first be imputed prior to fitting the STL and TBATS models, and this imputation scheme will have impacts on the resulting forecasts (though uncertainty in the imputation is somewhat respected); (2) the current modelling framework does not use partial pooling among locations and; (3) no environmental or ecological covariates are included. The purpose of the modelling framework is to identify useful benchmark models that can perform reasonably well in forecasting bounded ecological series such as tick densities so that more realistic population-based or mechanistic models that incorporate biological realism can be tested and improved.

## License

None whatsoever. Please use, adapt and distribute this code widely so that the models can be improved and our collective knowledge of tick abundance distributions can accelerate.
