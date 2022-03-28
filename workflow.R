#### Dependencies ####
# Data wrangling 
library(tidyverse)
library(MMWRweek)
library(lubridate)

# Modelling
library(forecast)
library(zoo)
library(runjags)
library(parallel)
library(MCMCvis)
library(mgcv)

# Metadata and forecast submission
library(EML)
library(neon4cast) #devtools::install_github('eco4cast/neon4cast')

#### Download tick data from EFI Github directly ####
read_csv("https://data.ecoforecast.org/targets/ticks/ticks-targets.csv.gz", 
         guess_max = 1e6) %>%
  dplyr::arrange(siteID, time) %>% 
  dplyr::mutate(year = lubridate::year(time)) -> targets

#### Fit Bayesian weighted ensemble benchmark models ####
source('R/benchmark_ens_jags.R')
unique_sites <- unique(targets$siteID)
last_obs_week <- 9
horizon = 4

# Plotting parameters
dir.create('Site_forecasts', recursive = T, showWarnings = F)
c_light <- c("#DCBCBC")
c_light_trans <- c("#DCBCBC70")
c_light_highlight <- c("#C79999")
c_mid <- c("#B97C7C")
c_mid_highlight <- c("#A25050")
c_mid_highlight_trans <- c("#A2505095")
c_dark <- c("#8F2727")
c_dark_highlight <- c("#7C0000")
probs <- c(0.05, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.95)

# Loop across sites and produce ensemble benchmark forecasts
site_fcs <- lapply(seq_along(unique_sites), function(site_index){
  cat('\n\nProcessing site', site_index, '...\n')
  
  # Create expanded site-specific dataframe with NAs for missing observations
  targets %>%
    dplyr::filter(siteID == unique_sites[site_index]) %>%
    dplyr::select(-time) %>%
    dplyr::arrange(year, mmwrWeek) -> targets_site
  targets_site %>%
    dplyr::full_join(expand.grid(mmwrWeek = seq(1, 52),
                                 year = seq(min(targets_site$year),
                                            2021),
                                 siteID = unique(targets_site$siteID))) %>%
    dplyr::arrange(year, mmwrWeek) -> targets_site
  targets_site[1:which(targets_site$mmwrWeek == last_obs_week &
                         targets_site$year == 2021), ] -> targets_site
  
  targets_site %>%
    dplyr::mutate(amblyomma_americanum = dplyr::case_when(
      mmwrWeek < 8 ~ 0,
      mmwrWeek > 45 ~ 0,
      TRUE ~ amblyomma_americanum
    )) -> targets_site
  
  # Convert outcome to ts object
  y <- ts(targets_site$amblyomma_americanum,
          frequency = 52, start = c(min(targets_site$year), 1))
  
  # Run the benchmark ensemble model to produce the forecast distribution
  site_fc <- benchmark_ens_jags(y, season = targets_site$mmwrWeek,
                                year = targets_site$year,
                                h = horizon)
  
  # Plot the forecast for a sanity check
  pdf(file = paste0('Site_forecasts/', unique_sites[site_index], '.pdf'),
      width = 6.5, height = 4.75)
  cred <- sapply(1:NCOL(site_fc),
                 function(n) quantile(site_fc[,n],
                                      probs = probs))
  plot(1, type = "n",
       xlab = 'Time',
       ylab = paste(unique_sites[site_index], 'forecast distribution'),
       xlim = c(1, length(y) + horizon),
       ylim = c(0, min(1000, max(cred))))
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
  points(as.vector(y), pch = 16, cex = 0.9, col = 'white')
  points(as.vector(y), pch = 16, cex = 0.65, col = 'black')
  abline(v = length(y), lty = 'dashed')
  dev.off()
  
  # Return 5000 samples from the posterior forecast distribution
  site_fc <- t(site_fc)
  t(site_fc[(length(y)+1):NROW(site_fc), sample(seq(1, NCOL(site_fc)),
                                                   5000, replace = F)])
  
})
names(site_fcs) <- unique_sites

# Save forecast distributions
dir.create('Forecasts', recursive = T, showWarnings = F)
save(site_fcs, file = 'Forecasts/site_fcs.rda')

#### Format forecast files ####
# All Sundays in 2021
sundays <- seq.Date(lubridate::ymd("2021-01-03"), 
                    by = 7, length.out = 52)

# Forecast mmwrWeeks
fx_weeks <- (last_obs_week + 1):(last_obs_week + 4)

# Forecast dates
fx_time <- sundays[fx_weeks]

# Gather forecasts into the correct format for the EFI Challenge Submission
load('Forecasts/site_fcs.rda')
all_submissions <- do.call(rbind, lapply(seq_along(unique_sites), function(x){
  site_fc_submission <- do.call(rbind, lapply(seq_len(NCOL(site_fcs[[x]])), function(fc_week){
    data.frame(time = MMWRweek::MMWRweek2Date(2021, fx_weeks[fc_week]),
               siteID = unique_sites[x],
               ensemble = seq(1, NROW(site_fcs[[x]])),
               forecast = 1,
               data_assimilation = 0,
               amblyomma_americanum = site_fcs[[x]][ , fc_week])
  }))
}))

# Save file as csv in the EFI format
# [theme_name]-[time]-[team_name].csv
theme_name <- "ticks"
time <- as.character(min(all_submissions$time))
team_name <- "TickBench"
file_name <- paste0(theme_name, "-", time, "-", team_name, ".csv.gz")
write_csv(all_submissions, file_name)

##### Format EML file ####
method_text <- paste(readLines('Metadata/methods.md'), 
                     collapse = "\n")
methods <- list(id = "forecast", 
                methodStep = list(description = 
                                    list(markdown = sub('epidemiological weeks', 
                                                        paste('epidemiological weeks weeks(', 
                                                              fx_weeks[1], '-', fx_weeks[4], 'in 2021)'), 
                                                        method_text))))
attributes <- tibble::tribble(
  ~attributeName,        ~attributeDefinition,               ~unit,           ~formatString,  ~numberType,    ~definition,
  "time",                 "time",                          "year",            "YYYY-MM-DD",  "numberType",     NA,
  "siteID",               "ID of NEON site",               "dimensionless",    NA,             NA,               "NEON site identifier",
  "ensemble",             "index of ensemble member",      "dimensionless",    NA,            "integer",          NA,
  "forecast",          "flag whether forecast or hindcast", "number",            NA,              "integer",       "flag whether a forecast or hindcast",
  "data_assimilation",    "flag whether includes assimilation", "dimensionless",NA, "integer",         NA,
  "amblyomma_americanum",    "forecast larval density of A. americanum","number",NA,    "real",           "forecast larval abundance of A. americanum"
)
attrList <- set_attributes(attributes, 
                           col_classes = c("Date", "character", 
                                           "numeric", "numeric", 
                                           "numeric", "numeric"))

physical <- set_physical(file_name)
dataTable <- eml$dataTable(
  entityName = file_name,
  entityDescription = "NEON tick abundance EFI forecast challenge - TickBench",
  physical = physical,
  attributeList = attrList)

me <- list(individualName = list(givenName = "Nicholas", 
                                 surName = "Clark"),
           electronicMailAddress = "n.clark@uq.edu.au",
           id = "https://orcid.org/0000-0001-7131-3301")

coverage <- 
  EML::set_coverage(begin = min(all_submissions$time), 
                    end = max(all_submissions$time),
                    sci_names = c('Amblyomma americanum'),
                    geographicDescription = "North American NEON tick survey sites",
                    westBoundingCoordinate = -124.848974,
                    eastBoundingCoordinate = 24.396308,
                    northBoundingCoordinate = 49.384358,
                    southBoundingCoordinate = -66.885444)

keywordSet <- list(
  list(
    keywordThesaurus = "EFI controlled vocabulary",
    keyword = list("forecast",
                   "amblyomma",
                   "ensemble",
                   "generalised additive model",
                   "bayesian computation",
                   "timeseries")))

dataset = eml$dataset(
  title = "NEON tick abundance EFI forecast challenge - TickBench",
  creator = me,
  contact = list(references="https://orcid.org/0000-0001-7131-3301"),
  pubDate = Sys.Date(),
  intellectualRights = "https://data.neonscience.org.",
  abstract =  "A Bayesian log-Gamma weighted ensemble forecast of NEON tick density.",
  dataTable = dataTable,
  keywordSet = keywordSet,
  coverage = coverage,
  methods = methods)

my_eml <- eml$eml(dataset = dataset,
                  packageId = "TickBench",  
                  system = "uuid")
EML::write_eml(my_eml, paste0(theme_name, "-", time, "-", team_name, ".xml"))

#### Submit using edited versions of the neon4cast scripts ####
source('R/submission_scripts.R')
submit(forecast_file = file_name, 
                  metadata = paste0(theme_name, "-", time, "-", team_name, ".xml"), 
                  ask = TRUE)
check_submission(file_name)

