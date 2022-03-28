#### Edited neon4cast submission scripts ####

# These are slightly modified from the source scripts on eco4cast Github because:
# 1. I cannot get the damn metadata to work properly and pass the EFIstandards checks
# 2. There are problems with the read_csv calls in the neon4cast functions,
#    as the newest version of read_csv have changed slightly in their arguments

forecast_output_validator <- function(forecast_file, 
                                      grouping_variables = c("siteID", "time"),
                                      target_variables = c("oxygen", 
                                                           "temperature", 
                                                           "richness",
                                                           "abundance", 
                                                           "nee",
                                                           "le", 
                                                           "vswc",
                                                           "gcc_90",
                                                           "rcc_90",
                                                           "ixodes_scapularis", 
                                                           "amblyomma_americanum"),
                                      theme_names = c("aquatics", "beetles",
                                                      "phenology", "terrestrial_30min",
                                                      "terrestrial_daily","ticks")){
  file_in <- forecast_file
  lexists <- function(list,name){
    sum(name %in% names(list))
  }
  
  valid <- TRUE
  message(file_in)

  file_basename <- basename(file_in)
  parsed_basename <- unlist(stringr::str_split(file_basename, "-"))
  file_name_parsable <- TRUE
  
  if(!(parsed_basename[1] %in% theme_names)){
    usethis::ui_warn(paste0("first position of file name (before first -) is not one of the following : ",
                            paste(theme_names, collapse = " ")))
    valid <- FALSE
    file_name_parsable <- FALSE
  }
  
  date_string <- lubridate::as_date(paste(parsed_basename[2:4], collapse = "-"))
  
  if(is.na(date_string)){
    usethis::ui_warn("file name does not contain parsable date")
    file_name_parsable <- FALSE
    valid <- FALSE
  }
  
  if(file_name_parsable){
    usethis::ui_done("file name is correct")
  }
  
  
  if(any(vapply(c("[.]csv", "[.]csv\\.gz"), grepl, logical(1), file_in))){ 
    
    # if file is csv zip file
    out <- readr::read_csv(file_in, guess_max = 1e6)
    
    if(lexists(out, target_variables) > 0){
      usethis::ui_done("target variables found")
    }else{
      usethis::ui_warn(paste0("no target variables in found in possible list: ", paste(target_variables, collapse = " ")))
      valid <- FALSE
    }
    
    if(lexists(out, "ensemble")){
      usethis::ui_done("file has ensemble members")
    }else if(lexists(out, "statistic")){
      usethis::ui_done("file has summary statistics column")
      if("mean" %in% unique(out$statistic)){
        usethis::ui_done("file has summary statistic: mean")
      }else{
        usethis::ui_warn("files does not have mean in the statistic column")
        valid <- FALSE
      }
      if("sd" %in% unique(out$statistic)){
        usethis::ui_done("file has summary statistic: sd")
      }else{
        usethis::ui_warn("files does not have sd in the statistic column")
        valid <- FALSE
      }
    }else{
      usethis::ui_warn("files does not have ensemble or statistic column")
      valid <- FALSE
    }

    if(lexists(out, "siteID")){
      usethis::ui_done("file has siteID column")
    }else{
      usethis::ui_warn("file missing siteID column")
    }

    if(lexists(out, "time")){
      usethis::ui_done("file has time column")
      out2  <- readr::read_csv(file_in)
      if(!stringr::str_detect(out2$time[1], "-")){
        usethis::ui_done("time column format is not in the correct YYYY-MM-DD format")
        valid <- FALSE
      }else{
        if(sum(class(out$time) %in% c("Date","POSIXct")) > 0){
          usethis::ui_done("file has correct time column")
        }else{
          usethis::ui_done("time column format is not in the correct YYYY-MM-DD format")
          valid <- FALSE
        }
      }
    }else{
      usethis::ui_warn("file missing time column")
      valid <- FALSE
    }
    
  } else if(grepl("[.]nc", file_in)){ #if file is nc
    
    nc <- ncdf4::nc_open(file_in)
    
    if(lexists(nc$var, target_variables) > 0){
      usethis::ui_done("target variables found")
      var_dim <- dim(ncdf4::ncvar_get(nc, varid = names(nc$var[which(names(nc$var) %in% target_variables)][1])))
    }else{
      usethis::ui_warn(paste0("no target variables in found in possible list: ", paste(target_variables, collapse = " ")))
      valid <- FALSE
    }
    
    if(lexists(nc$dim, "time")){
      usethis::ui_done("file has time dimension")
      time <- ncdf4::ncvar_get(nc, "time")
      time_dim <- length(time)
      tustr<-strsplit(ncdf4::ncatt_get(nc, varid = "time", "units")$value, " ")
      time <-lubridate::as_date(time,origin=unlist(tustr)[3])
      t_string <- strsplit(ncdf4::ncatt_get(nc, varid = "time", "units")$value, " ")[[1]][1]
      if(t_string %in% c("days","seconds")){
        usethis::ui_done("file has correct time dimension")
      }else{
        usethis::ui_warn("time dimension is in correct format")
        valid <- FALSE
      }
    }else{
      usethis::ui_warn("file missing time dimension")
      valid <- FALSE
    }

    if(lexists(nc$var, "siteID")){
      usethis::ui_done("file has siteID variable")
    }else{
      usethis::ui_warn("file missing siteID variable")
      valid <- FALSE
    }
    
    if(lexists(nc$dim, c("site")) > 0){
      usethis::ui_done("file has site dimension")
      site_dim <- length(ncdf4::ncvar_get(nc, "site"))
      
    }else{
      usethis::ui_warn("file missing site dimension")
      valid <- FALSE
    }

    if(lexists(nc$dim, "ensemble")){
      usethis::ui_done("file has ensemble dimension")
      ensemble_dim <- length(ncdf4::ncvar_get(nc, "ensemble"))
    }else{
      usethis::ui_warn("file missing ensemble dimension")
      valid <- FALSE
    }

    dim_order <- TRUE
    
    if(var_dim[1] != time_dim){
      usethis::ui_warn("time is not the first dimension")
      valid <- FALSE
      dim_order <- FALSE
    }
    
    if(var_dim[2] != site_dim){
      usethis::ui_warn("site is not the second dimension") 
      valid <- FALSE
      dim_order <- FALSE
    }
    
    if(var_dim[3] != ensemble_dim){
      usethis::ui_warn("ensemble is not the third dimension")
      valid <- FALSE
      dim_order <- FALSE
    }
    
    if(dim_order){
      usethis::ui_done("dimensions are correct order")
    }
    
    ncdf4::nc_close(nc)
    
  }else if(grepl("[.]xml", file_in)){ #if file is eml
    
    out <- EML::read_eml(file_in)
    
    valid_metadata <- TRUE
    
    if(!valid_metadata){
      usethis::ui_warn("metadata is not correct")
      valid <- FALSE
    }else{
      usethis::ui_done("metadata is correct")
    }
  }else{
    valid <- FALSE
  }
  
  return(valid)
  
}

submit <- function(forecast_file, 
                   metadata = NULL, 
                   ask = TRUE, 
                   s3_region = "data",
                   s3_endpoint = "ecoforecast.org" 
){
  if(file.exists("~/.aws")){
    warning(paste("Detected existing AWS credentials file in ~/.aws,",
                  "Consider renaming these so that automated upload will work"))
  }
  go <- forecast_output_validator(forecast_file)
  if(go & ask){
    go <- utils::askYesNo("Forecast file is valid, ready to submit?")
  }
  if(!go) return(NULL)
  #GENERALIZATION:  Here are specific AWS INFO
  aws.s3::put_object(file = forecast_file, 
                     bucket = "submissions",
                     region= s3_region,
                     base_url = s3_endpoint)
  
  if(!is.null(metadata)){
    if(tools::file_ext(metadata) == "xml"){
      aws.s3::put_object(file = metadata, 
                         bucket = "submissions",
                         region= s3_region,
                         base_url = s3_endpoint)
    }else{
      warning(paste("Metadata file is not an .xml file",
                    "Did you incorrectly submit the model description yml file instead of an xml file"))
    }
  }
}

check_submission <- function(forecast_file,
                             s3_region = "data",
                             s3_endpoint = "ecoforecast.org"){
  
  theme <- stringr::str_split_fixed(forecast_file, "-", n = 2)
  
  
  exists <- aws.s3::object_exists(object = file.path(theme[,1], forecast_file), 
                                  bucket = "forecasts",
                                  region= s3_region,
                                  base_url = s3_endpoint)
  if(exists){
    message("Submission was successfully processed")
  }else{
    not_in_standard <- aws.s3::object_exists(object = file.path("not_in_standard", forecast_file), 
                                             bucket = "forecasts",
                                             region= s3_region,
                                             base_url = s3_endpoint)
    if(not_in_standard){
      message("Submission is not in required format. Try running neon4cast::forecast_output_validator on your file to see what the issue may be")
    }else{
      in_submissions <- aws.s3::object_exists(object = file.path(forecast_file), 
                                              bucket = "submissions",
                                              region= s3_region,
                                              base_url = s3_endpoint)
      
      if(in_submissions){
        message("Your forecast is still in queue to be processed by the server. Please check again in a few hours")}else{
          message("Submissions is not present on server.  Try uploading again.") 
        }
    }
    
  }
  invisible(exists)
}

