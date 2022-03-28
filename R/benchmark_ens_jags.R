benchmark_ens_jags = function(y, season, year, h = 4){
  
  # Impute any NAs in y using a robust seasonal imputation
  y_impute <- as.vector(forecast::na.interp(y))
  y_impute <- ts(y_impute, frequency = 52, start = c(min(targets_site$year), 1))
  y_impute[y_impute < 0] <- 0
  
  # Log(y + 1/6) for the two benchmark timeseries models that require Gaussian responses
  # (STLM and TBATS)
  log_y_impute <- log(y_impute + 1/6)
  
  # Drop the last (h + 26) observations so that each benchmark's forecast can be somewhat evaluated
  # when quantifying ensemble weights
  indices_train <- 1:(length(y) - (h + 26))
  log_y_train <- ts(log_y_impute[indices_train], start = start(y_impute),
                    frequency = frequency(y_impute))
  
  #### Fit benchmark time series models to log(x + 1/6) training series and forecast
  # the full out-of-sample period ####
  # STLM with AR residuals
  fc_mod <- forecast::stlm(log_y_train, modelfunction = ar)
  fc_fc <- forecast::forecast(fc_mod, (h * 2) + 26)
  
  # Extract fitted in-sample and forecasted out-of-sample values
  stlm_fit <- as.vector(forecast::na.interp(zoo::na.approx(c(fc_mod$fitted, 
                                                             fc_fc$mean), na.rm = F)))
  
  # Repeat for TBATS
  fc_mod <- forecast::tbats(log_y_train)
  fc_fc <- forecast::forecast(fc_mod, (h * 2) + 26)
  tbats_fit <- c(fc_mod$fitted.values, fc_fc$mean)
  
  # Repeat for a GAM model that uses smooths to jointly model the location and scale of a
  # Gamma distribution for the raw in-sample data
  gam_mod <- gam(list(y ~ s(season, bs = 'cc', k = 12) + s(year, bs = 'gp', k = 4), 
                        ~ s(season, bs = 'cc', k = 12) + s(year, bs = 'gp', k = 4)),
                      family = gammals,
                 data = data.frame(y = (y + 0.0001)[indices_train], 
                                   season = season[indices_train], 
                                   year = year[indices_train]))
  gam_fit <- predict(gam_mod, type = 'response', newdata = rbind(data.frame(season = season, year = year),
                                                      data.frame(season = seq(tail(season, 1),
                                                                              tail(season, 1) + (h-1)),
                                                                 year = 2021)))
  gam_fit <- log(rgamma(n = NROW(gam_fit), 
                        shape = gam_fit[,1], 
                        scale = exp(gam_fit[,2])) + 1/6)
  
  #### Use fitted values from each benchmark as weighted predictors for forecasting the raw
  # original, unimputed time series ####
  y_orig <- y
  y <- c(as.vector(y), rep(NA, h)) + 0.0001
  
  model_code <- " model {
  
  ## Mean expectation (shape) linear predictor
  for (i in 1:n) {
   mu[i] <- exp(alpha + 
                beta[1] * stlm_fit[i] + 
                beta[2] * gam_fit[i] +
                beta[3] * tbats_fit[i] +
                trend[i])
  }

  ## Random walk trend to capture any remaining temporal autocorrelation
  trend[1] ~ dnorm(0, tau)
  for (i in 2:n){
   trend[i] ~ dnorm(trend[i - 1], tau)
  }

  ## Gamma likelihood
  for (i in 1:n) {
   y[i] ~ dgamma(mu[i], lambda)
  }
  
  ## Posterior predictions
  for (i in 1:n) {
   ypred_raw[i] ~ dgamma(mu[i], lambda);
   ypred[i] <- max(ypred_raw[i] - 0.0001, 0)
  }
  
  ## Priors, accounting for likely correlation among model weights
  beta[1:3] ~ dmnorm(zeros, Omega)
  Omega ~ dscaled.wishart(scales, 2)
  zeros <- c(0, 0, 0)
  alpha ~ dunif(-1, 1)
  lambda ~ dexp(1)
  sigma ~ dexp(2.5)
  tau <- pow(sigma, -2)
}"

  jags_data <- list(
    y = y,
    scales = rep(1, 3),
    stlm_fit = stlm_fit,
    gam_fit = gam_fit,
    tbats_fit = tbats_fit,
    n = length(y))
  
  # Run the model using 4 parallel MCMC chains in runjags
  chains <- 4
  cl <- parallel::makePSOCKcluster(min(c(chains, parallel::detectCores() - 1)))
  setDefaultCluster(cl)
  unlink('jags_mod.txt')
  cat(model_code, file = 'jags_mod.txt', sep = '\n', append = T)
  model_run <- run.jags(
    data = jags_data,
    n.chains = chains,
    modules = 'glm',
    adapt = 10000,
    burnin = 5000,
    sample = 5000,
    thin = 5,
    method = "rjparallel",
    monitor = c('ypred',
                'sigma',
                'beta'),
    cl = cl,
    model = 'jags_mod.txt'
  )
  stopCluster(cl)
  model_run <- coda::as.mcmc.list(model_run)
  unlink('jags_mod.txt')
  
  # Visualise weight posteriors as a very simple sanity check of convergence
  MCMCvis::MCMCtrace(model_run, c('beta'), pdf = F, n.eff = T,
                     main_den = c('stlm', 'gam', 'tbats'),
                     main_tr = c('stlm', 'gam', 'tbats'))
  
  # Return forecast distribution
  ypreds <- MCMCvis::MCMCchains(model_run, 'ypred')
  ypreds
}

#### Function to plot posterior distribution ####
plot_posterior = function(site_fc, y, horizon, 
                          sitename = unique_sites[site_index],
                          year){
  
  # Colour scheme
  c_light <- c("#DCBCBC")
  c_light_trans <- c("#DCBCBC70")
  c_light_highlight <- c("#C79999")
  c_mid <- c("#B97C7C")
  c_mid_highlight <- c("#A25050")
  c_mid_highlight_trans <- c("#A2505095")
  c_dark <- c("#8F2727")
  c_dark_highlight <- c("#7C0000")
  
  # Posterior empirical quantiles
  probs <- c(0.05, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.95)
  cred <- sapply(1:NCOL(site_fc),
                 function(n) quantile(site_fc[,n],
                                      probs = probs))
  
  # Plot
  plot(1, type = "n",
       ylab = paste(sitename, 'posterior predictive distribution'),
       xlim = c(1, length(y) + horizon),
       ylim = c(0, min(1000, max(cred))),
       xlab = '',
       xaxt = 'n')
  pred_vals <- seq(1, length(y) + horizon) 
  polygon(c(pred_vals, rev(pred_vals)), c(cred[1,], rev(cred[9,])),
          col = c_light, border = NA)
  polygon(c(pred_vals, rev(pred_vals)), c(cred[2,], rev(cred[8,])),
          col = c_light_highlight, border = NA)
  polygon(c(pred_vals, rev(pred_vals)), c(cred[3,], rev(cred[7,])),
          col = c_mid, border = NA)
  polygon(c(pred_vals, rev(pred_vals)), c(cred[4,], rev(cred[6,])),
          col = c_mid_highlight, border = NA)
  lines(pred_vals, cred[5,], col = c_dark, lwd = 2.5)
  points(as.vector(y), pch = 16, cex = 0.75, col = 'white')
  points(as.vector(y), pch = 16, cex = 0.5, col = 'black')
  abline(v = length(y), lty = 'dashed')
  axis(1, at = seq(0, length(y) + horizon,
                   b = 52), labels = seq(min(year), max(year)), cex.axis = 1)
}


