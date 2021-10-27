# clear workspace
rm(list = ls())

# packages
packs <- c("dplyr", "nimble", "htmltools", "ggplot2", "sf", "Rcpp", "RcppArmadillo", "inline", "mvtnorm", "readr", "parallel", "xtable", "rstan", "coda", "vegan", "tidyverse", "lubridate", "bcmaps", "knitr", "kableExtra", "ggalluvial")
sapply(packs, require, character.only = T)
rm(packs)

# convenience
`%notin%` <- Negate("%in%")
options(mc.cores = parallel::detectCores())

######################
### Clean response ###
######################
get_night <- function(ts){
  year <- lubridate::year(ts)
  month <- lubridate::month(ts)
  day <- lubridate::day(ts)
  hour <- lubridate::hour(ts)
  
  tmp <- tibble(
    year = year,
    month = month,
    day = day,
    hour = hour
  ) %>%
    mutate(
      date = ymd_h(paste(paste(year, month, day, sep = "-"), hour))
    ) %>%
    mutate(
      sample_night = case_when(
        hour >= 12 & hour <= 24 ~ ymd(paste(year, month, day, sep = "-")),
        TRUE ~ ymd(paste(year, month, day, sep = "-")) - 1
      )
    )
  
  return(tmp$sample_night)
}

# clean raw data
bcdata <- readRDS("bcdata.rds")

# remove files with no auto ID and transect files
bcdata2 <- bcdata %>%
  filter(`WA|Kaleidoscope|Auto ID` %notin% c("NoID", "NOISE")) %>%
  filter(Quadrant != "transects", Quadrant != "Transects")

# figure out sample night - replace year in Timestamp with Year field as per Jason's email
bcdata2 <- bcdata2 %>%
  mutate(
    year_ts = lubridate::year(Timestamp),
    month_ts = lubridate::month(Timestamp),
    day_ts = lubridate::day(Timestamp),
    hour_ts = lubridate::hour(Timestamp),
    min_ts = lubridate::minute(Timestamp),
    sec_ts = lubridate::second(Timestamp),
    ts = lubridate::ymd_hms(paste(paste(Year, month_ts, day_ts, sep = "-"), paste(hour_ts, min_ts, sec_ts, sep = "-")))
  ) %>%
  select(ts, everything()) %>%
  mutate(
    sample_night = get_night(ts)
  ) %>%
  select(ts, Timestamp, sample_night, everything()) %>%
  select(-ends_with("ts"))

bcclean <- bcdata2 %>%
  mutate(
    man = case_when(
      nchar(`Final Species ID`) %notin% c(4, 5) ~ "other",
      grepl("Q", `Final Species ID`) ~ "other",
      grepl("[[:digit:]]", `Final Species ID`) ~ "other",
      grepl("High", `Final Species ID`) ~ "other",
      grepl("HIGH", `Final Species ID`) ~ "other",
      grepl("Hi", `Final Species ID`) ~ "other",
      grepl("LOW", `Final Species ID`) ~ "other",
      grepl("Low", `Final Species ID`) ~ "other",
      grepl("Bat", `Final Species ID`) ~ "other",
      grepl("BAT", `Final Species ID`) ~ "other",
      TRUE ~ `Final Species ID`
    )
  ) %>%
  mutate(man = case_when(
    substring(man, 1, 1) %in% c("k", "m") ~ substring(man, 2), 
    TRUE ~ man
  )) %>%
  select(
    sample_night = sample_night,
    year = Year,
    grts = `GRTS Cell ID`,
    quad = Quadrant,
    man_dirty = `Final Species ID`,
    man = man, 
    kal = `WA|Kaleidoscope|Auto ID`,
    sb = `SB|Species Auto ID`
  )

#######################
### Site covariates ###
#######################
# load in .shp from ssa
site_covariates <- readRDS("site_covariates.rds") %>%
  st_as_sf()

# get BC boundary
bc_boundary <- bcmaps::bc_bound("sf") %>%
  st_transform(crs = 4269)

# trim ssa files to bc_boundary
bc_nabat <- site_covariates %>%
  slice(
    st_intersects(bc_boundary, .) %>% unlist() %>% unique() %>% sort
  ) %>%
  select(
    grts = grts_id, 
    mean_elevation,
    annual_precip,
    annual_mean_temp,
    percent_forest,
    percent_water,
    geometry
  ) 

# join
bc_join <- left_join(
  bcclean, 
  bc_nabat, 
  by = "grts"
) %>%
  filter(man %notin% c("other", "NOISE", "NoID", "Noise", "noise", "noID")) %>%
  filter_at(vars(mean_elevation:percent_water), all_vars(!is.na(.)))

###########################
### Activity covariates ###
###########################
# get locations and times
daymet <- bc_join %>%
  dplyr::select(sample_night, year, grts, geometry) %>%
  distinct %>%
  st_as_sf() %>%
  mutate(centroid = st_geometry(st_centroid(.))) %>%
  as_tibble %>%
  dplyr::select(-geometry) %>%
  rename(geometry = centroid) %>%
  st_as_sf() %>%
  mutate(
    lon = st_coordinates(.)[,1],
    lat = st_coordinates(.)[,2]
  ) %>%
  as_tibble

daymet_nabat <- function(tbl){
  
  out <- list()
  for(i in 1:nrow(tbl)){
    
    tryCatch({
      message(paste0("Scraping DayMet data for row ", i))
      tmp <- suppressMessages(daymetr::download_daymet(site = tbl$grts[i],
                                                       lat = tbl$lat[i],
                                                       lon = tbl$lon[i],
                                                       start = tbl$year[i],
                                                       end = tbl$year[i])$data)
      
      # rename columns
      names(tmp) <- c("year", "year_day", "day_length", "precip", "shortwave_radiation", 
                      "snow_water_equiv", "max_air_temp", "min_air_temp", "water_vap_pressure")
      
      # add GRTS cell identifier
      tmp$grts <- rep(tbl$grts[i], nrow(tmp))
      
      # get year-day
      year_day_target <- lubridate::yday(tbl$sample_night[i])
      
      tmp <- filter(tmp, year_day == year_day_target) %>%
        mutate(month = month(tbl$sample_night[i]), day = day(tbl$sample_night[i]), sample_night = tbl$sample_night[i])
      
      # store
      out[[i]] <- tmp
      
      if(i %% 100 == 0){
        saveRDS(out, file = "processed covariates/out.rds")
      }
      
      message(paste0("Progress: ", i, "/", nrow(tbl)))
    }, error = function(e) {out[[i]] <- paste0(conditionMessage(e), "\n")})
    
  }
  
  return(out)
  
}
daymet_out <- daymet_nabat(daymet)
saveRDS(daymet_out, file = "daymet_out.rds")

# find missing locations - they're all in water!
missing_daymet <- daymet[which(sapply(daymet_out, is.null)),] 
missing_daymet_sf <- st_as_sf(missing_daymet, coords = c("lon", "lat"), crs = 4269)

# select point inside cell till it works
daymet_out2 <- list()
for(i in 1:nrow(missing_daymet)){
  # isolate row
  tmp <- missing_daymet[i,]
  cell <- bc_join[which(bc_join$grts == tmp$grts), ] %>% select(grts, geometry) %>% distinct() %>% st_as_sf()
  
  signal <- FALSE
  while(!signal){
    rand_point <- suppressMessages(st_coordinates(st_sample(cell, size = 1)))
    tmp2 <- tryCatch({
      suppressMessages({
        daymetr::download_daymet(
          lat = rand_point[2], lon = rand_point[1], start = tmp$year, end = tmp$year
        )
      })
    }, error = function(e){FALSE})
    
    if(class(tmp2) == "daymetr"){signal <- TRUE}
  }
  
  # tmp2 <- tmp2$data
  names(tmp2$data) <- c("year", "year_day", "day_length", "precip", "shortwave_radiation",
                        "snow_water_equiv", "max_air_temp", "min_air_temp", "water_vap_pressure")
  
  # post process
  daymet_out2[[i]] <- tmp2$data %>%
    filter(year_day == lubridate::yday(tmp$sample_night)) %>%
    mutate(
      month = month(tmp$sample_night),
      day = day(tmp$sample_night),
      sample_night = tmp$sample_night,
      grts = tmp$grts
    ) %>%
    dplyr::select(grts, year, month, day, everything())
  
  # progress
  message(paste0("Progress: ", i, "/", nrow(missing_daymet)))
}
saveRDS(daymet_out2, "daymet_out2.rds")

# bind list into dataframe
bc_final <- bind_rows(
  do.call("bind_rows", daymet_out),
  do.call("bind_rows", daymet_out2)
) %>%
  dplyr::select(grts, sample_night, year, month, day, everything()) %>%
  left_join(dplyr::select(bc_join, -year), ., by = c("grts", "sample_night"))


# fix year
bc_final <- bc_final %>%
  mutate(year = year(sample_night))

# add lunar illumination
bc_final <- bc_final %>%
  mutate(
    lunar_illum = lunar::lunar.illumination(
      x = sample_night
    )
  )

# export
saveRDS(bc_final, "bc_final.rds")