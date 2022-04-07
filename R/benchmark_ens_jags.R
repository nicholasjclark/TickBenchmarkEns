benchmark_ens_jags = function(y, season, year, h = 4, num_draws = 10, chains = 4){
  
  # Set cross-validation indices
  indices_train <- 1:(length(y) - (h + 52))
  
  # Fit a GAM with Tweedie distributed response
  gam_mod <- gam(y ~ s(season, bs = 'cc', k = 8) + 
                   s(year, bs = 'bs', k = 5, m = c(2, 1)) +
                   ti(season, year, bs = c('cc', 'gp'), k = c(8, 3)),
                 knots = list(season = c(0.5, 52.5), 
                              year = c(min(year) - 1, 
                                       min(year), 
                                       max(year), 
                                       max(year) + 1)),
                 family = tw,
                 data = data.frame(y = (y)[indices_train], 
                                   season = season[indices_train], 
                                   year = year[indices_train]))
  
  # Extract posterior GAM expectations for training and testing periods
  rmvn <- function(n,mu,sig){
    L <- mroot(sig); m <- ncol(L);
    t(mu + L %*% matrix(rnorm(m*n), m, n)) 
  }
  betas <- rmvn(num_draws, coef(gam_mod), gam_mod$Vp) 
  Xp <- predict(gam_mod, newdata = data.frame(season = c(season, 
                                                         seq(tail(season, 1) + 1,
                                                           tail(season, 1) + (h))), 
                                              year = c(year,
                                                       rep(2021, h))), 
                type = "lpmatrix") 
  gam_fits <- matrix(NA, nrow = num_draws, ncol = (length(y)+h))
  for(i in 1:num_draws){ 
    gam_fits[i, ] <- Xp %*% betas[i,]
  }
  
  # Extract sarima of GAM-imputed training values
  gam_impute <- function(mu, y){
    # If y is NA, use a different posterior expectation
    y_working <- vector()
    for(i in 1:(length(y)+h)){
      y_working[i] <- ifelse(is.na(y[i]), mu[i], log(y[i] + 0.0001))
    }
    y_working 
  }
  
  gam_imputed <- matrix(NA, nrow = num_draws, ncol = (length(y) + h))
  for(i in 1:num_draws){
    gam_imputed[i, ] <- gam_impute(mu = gam_fits[i,],
                                   y = y)
  }

  #### Fit benchmark time series models to log-imputed training series and forecast
  # the full out-of-sample period ####
  arima_mod <- forecast::auto.arima(ts(gam_imputed[1, ], 
                                       start = c(min(year), 1),
                                       frequency = 52),
                                    nmodels = 40, approximation = TRUE)
  cl <- parallel::makePSOCKcluster(min(c(chains, parallel::detectCores() - 2)))
  setDefaultCluster(cl)
  clusterExport(NULL, c('gam_imputed',
                        'year',
                        'h',
                        'y',
                        'indices_train',
                        'arima_mod'),
                envir = environment())
  clusterEvalQ(cl, library(forecast))
  clusterEvalQ(cl, library(zoo))
  
  pbapply::pboptions(type = "none")
  fc_mod_fits <- pbapply::pblapply(seq_len(num_draws), function(i){
    imputed_ts <- ts(gam_imputed[i, indices_train], 
                       start = c(min(year), 1),
                       frequency = 52)
    
    # STLM with AR residuals
    fc_mod <- forecast::stlm(imputed_ts, modelfunction = ar)
    fc_fc <- forecast::forecast(fc_mod, (h * 2) + 52)
    # Extract fitted in-sample and forecasted out-of-sample values
    stlm_fits <- as.vector(forecast::na.interp(zoo::na.approx(c(fc_mod$fitted, 
                                                                         fc_fc$mean), na.rm = F)))
    
    # Repeat for SARIMA
    fc_mod <- forecast::Arima(imputed_ts, model = arima_mod)
    fc_fc <- forecast::forecast(fc_mod, (h * 2) + 52)
    arima_fits <- c(fc_mod$fitted, fc_fc$mean)
    list(stlm_fits = stlm_fits,
         arima_fits = arima_fits)
    
  }, cl = cl)
  stopCluster(cl)
  
  stlm_fits <- do.call(rbind, purrr::map(fc_mod_fits, 'stlm_fits'))
  arima_fits <- do.call(rbind, purrr::map(fc_mod_fits, 'arima_fits'))
  rm(fc_mod_fits)
  
  #### Use fitted values from each benchmark as weighted predictors for forecasting the raw
  # original time series ####
  # Add a small offset so there are no exact zeros; won't impact forecasts
  # or inference but will stop problems due to starting values for the Gamma
  # parameters
  y_orig <- y
  y <- c(as.vector(y), rep(NA, h)) + 0.0001
  
  # Forecast values get higher weights than fitted (training) values
  weight <- c(rep(1, length(indices_train)),
              rep(1.5, (h * 2) + 52))
  
  model_code <- " model {
  
  ## Mean expectation (shape) linear predictor
  for (i in 1:n) {
   mu[i] <- exp(beta[1] * gam_fits[K, i] + 
                beta[2] * stlm_fits[K, i] + 
                beta[3] * arima_fits[K, i])
  }
  
  ## Likelihood function
  for (i in 1:n) {
   y[i] ~ dgamma(mu[i], lambda * weight[i])
  }
  
  ## Posterior predictions
  for (i in 1:n) {
   ypred[i] ~ dgamma(mu[i], lambda * weight[i])
  }
  
  ## Priors
  # Account for likely correlation among model weights
  lambda ~ dexp(3)
  beta[1:3] ~ dmnorm(zeros, Omega)
  # Rename betas for easier extraction
  beta_gam <- beta[1]
  beta_stlm <- beta[2]
  beta_arima <- beta[3]
  Omega ~ dscaled.wishart(scales, 2)
  zeros <- c(0, 0, 0)
  
  # Design matrices (sample with equal probabilities)
  K ~ dcat(design_probs)
  for (k in 1:num_draws) {
    design_probs[k] <- 1
  }
}"

  jags_data <- list(
    y = y,
    scales = rep(1, 3),
    stlm_fits = stlm_fits,
    gam_fits = gam_fits,
    arima_fits = arima_fits,
    num_draws = num_draws,
    weight = weight,
    n = length(y))
  
  # Run the model using 4 parallel MCMC chains in runjags
  initlist <- replicate(chains, list(beta = runif(3, -1, 1)),
                        simplify = FALSE)
  cl <- parallel::makePSOCKcluster(min(c(chains, parallel::detectCores() - 1)))
  setDefaultCluster(cl)
  unlink('jags_mod.txt')
  cat(model_code, file = 'jags_mod.txt', sep = '\n', append = T)
  model_run <- run.jags(
    data = jags_data,
    n.chains = chains,
    modules = 'glm',
    adapt = 2000,
    burnin = 2000,
    sample = 2000,
    thin = 2,
    inits = initlist,
    method = "rjparallel",
    monitor = c('ypred',
                'lambda',
                'beta_gam',
                'beta_stlm',
                'beta_arima'),
    cl = cl,
    model = 'jags_mod.txt')
  stopCluster(cl)
  model_run <- coda::as.mcmc.list(model_run)
  unlink('jags_mod.txt')
  
  # Visualise weight posteriors as a very simple sanity check of convergence
  MCMCvis::MCMCtrace(model_run, c('beta_gam',
                                  'beta_stlm',
                                  'beta_arima',
                                  'lambda'), pdf = F, n.eff = T)
  
  #### Refit all benchmark models to the entire observation history and produce forecasts;
  # then weight forecasts according to the estimated weights from the JAGS model ####
  # Re-estimated GAM
  gam_mod <- gam(y ~ s(season, bs = 'cc', k = 8) + 
                   s(year, bs = 'bs', k = 5, m = c(2, 1)) +
                   ti(season, year, bs = c('cc', 'gp'), k = c(8, 3)),
                 knots = list(season = c(0.5, 52.5), 
                              year = c(min(year) - 1, 
                                       min(year), 
                                       max(year), 
                                       max(year) + 1)),
                 family = tw,
                 data = data.frame(y = y_orig, 
                                   season = season, 
                                   year = year))
  betas <- rmvn(num_draws, coef(gam_mod), gam_mod$Vp) 
  Xp <- predict(gam_mod, newdata = data.frame(season = c(season, 
                                                         seq(tail(season, 1) + 1,
                                                             tail(season, 1) + (h))), 
                                              year = c(year,
                                                       rep(2021, h))), 
                type = "lpmatrix") 
  gam_fits <- matrix(NA, nrow = num_draws, ncol = (length(y_orig)+h))
  for(i in 1:num_draws){ 
    gam_fits[i, ] <- Xp %*% betas[i,]
  }
  gam_imputed <- matrix(NA, nrow = num_draws, ncol = (length(y_orig) + h))
  for(i in 1:num_draws){
    gam_imputed[i, ] <- gam_impute(mu = gam_fits[i,],
                                   y = y_orig)
  }
  
  cl <- parallel::makePSOCKcluster(min(c(chains, parallel::detectCores() - 2)))
  setDefaultCluster(cl)
  clusterExport(NULL, c('gam_imputed',
                        'year',
                        'h',
                        'y_orig',
                        'arima_mod'),
                envir = environment())
  clusterEvalQ(cl, library(forecast))
  clusterEvalQ(cl, library(zoo))
  
  pbapply::pboptions(type = "none")
  fc_mod_fits <- pbapply::pblapply(seq_len(num_draws), function(i){
    imputed_ts <- ts(gam_imputed[i, 1:length(y_orig)], 
                       start = c(min(year), 1),
                       frequency = 52)
    
    # STLM with AR residuals
    fc_mod <- forecast::stlm(imputed_ts, modelfunction = ar)
    fc_fc <- forecast::forecast(fc_mod, h)
    # Extract fitted in-sample and forecasted out-of-sample values
    stlm_fits <- as.vector(forecast::na.interp(zoo::na.approx(c(fc_mod$fitted, 
                                                                     fc_fc$mean), na.rm = F)))
    
    # Repeat for arima
    fc_mod <- forecast::Arima(imputed_ts, model = arima_mod)
    fc_fc <- forecast::forecast(fc_mod, h)
    arima_fits <- c(fc_mod$fitted, fc_fc$mean)
    list(stlm_fits = stlm_fits,
         arima_fits = arima_fits)
    
  }, cl = cl)
  stopCluster(cl)
  
  stlm_fits <- do.call(rbind, purrr::map(fc_mod_fits, 'stlm_fits'))
  arima_fits <- do.call(rbind, purrr::map(fc_mod_fits, 'arima_fits'))
  rm(fc_mod_fits)
  
  # Extract posterior parameter estimates
  beta_gams <- MCMCvis::MCMCchains(model_run, 'beta_gam')
  beta_stlms <- MCMCvis::MCMCchains(model_run, 'beta_stlm')
  beta_arimas <- MCMCvis::MCMCchains(model_run, 'beta_arima')
  lambdas <- MCMCvis::MCMCchains(model_run, 'lambda')
  
  ypreds <- do.call(rbind, lapply(seq_len(2000), function(x){
    # Sample indices
    design_index <- sample(1:num_draws, 1, T)
    param_index <- sample(1:NROW(beta_gams), 1, T)
    
    # Sampled design matrix
    gam_fits_samp <- gam_fits[design_index, ]
    stlm_fits_samp <- stlm_fits[design_index, ]
    arima_fits_samp <- arima_fits[design_index, ]
    
    # Sampled parameters
    b_gam <- beta_gams[param_index,]
    b_stlm <- beta_stlms[param_index,]
    b_arima <- beta_arimas[param_index,]
    lambda <- lambdas[param_index,]
    
    # Posterior expectations
    mu <- exp(b_gam * gam_fits_samp + 
              b_stlm * stlm_fits_samp + 
              b_arima * arima_fits_samp)
    
    # Posterior predictions (include the weights for forecasted values)
    rgamma(length(mu),
           shape = mu, 
           rate = lambda * 1.5)
  }))
  
  # Return forecast distribution
  ypreds
}

#### Function to plot posterior distribution ####
plot_posterior = function(site_fc, y, horizon, 
                          sitename = unique_sites[site_index],
                          year){
  
  .pardefault <- par(no.readonly=T)
  par(.pardefault)
  par(mfrow = c(1, 1))
  
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
  par(.pardefault)
}


