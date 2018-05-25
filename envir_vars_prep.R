#Clean and prep meteorology, pollution variables
rm(list=ls())
#set up smart sensing of os
if(Sys.info()[1]=='Windows'){
  j = "J:/"
} else{
  j = '/home/j/'
}

# load packages, install if missing  
pack_lib = '/snfs2/HOME/wgodwin/R'
.libPaths(pack_lib)
pacman::p_load(data.table, fst, ggplot2, parallel, magrittr, lavaan, foreign)

#set settings
main.dir <- paste0(j, "temp/wgodwin/thesis/data/01_raw/")
save.dir <- paste0(j, "temp/wgodwin/thesis/data/02_prepped/")
model.dir <- paste0(j, "temp/wgodwin/thesis/data/03_analysis/")
meta.dir <- paste0(j, "temp/wgodwin/thesis/metadata/")
pop.dir <- paste0(j, "WORK/05_risk/risks/temperature/data/exp/pop/mex/all_age/")
muni.dt <- fread(paste0(meta.dir, "station_muni_map.csv"))
day <- F
######################################Function to calculate moving average over 30 day lag############################
mean_lag_3 <- function(i, data, var, lag) {
  is.near <- as.numeric(data$date[i] - data$date) >= 0 & as.numeric(data$date[i] - data$date) < lag
  #mean(data$paste(var)[is.near], na.rm = T)
  data[, mean(get(var)[is.near], na.rm = T)]
}

##################################################################
#######temperature prep
##################################################################
dt.all <- data.table()
for(year in 1998:2016){
  #read in data and subset to temperature parameter
  dt <- fread(paste0(main.dir, "meteorology/meteorología_", year, ".csv"))
  setnames(dt, c(colnames(dt)[grep("station", colnames(dt))], colnames(dt)[grep("parameter", colnames(dt))], "value"),
           c("id_station", "parameter", "temp_mean_day")) #min, max, mean
  dt2 <- fread(paste0(main.dir, "meteorology/meteorología_", paste(year -1), ".csv"))
  setnames(dt2, c(colnames(dt2)[grep("station", colnames(dt2))], colnames(dt2)[grep("parameter", colnames(dt2))], "value"),
           c("id_station", "parameter", "temp_mean_day")) #min, max, mean
  dt <- rbind(dt, dt2)
  rm(dt2)
  dt <- dt[parameter == "TMP"] #RH, TMP, WSP
  
  #Read in year previous
  first <- as.Date(paste0("01/", as.numeric(year)), format = "%j/%Y")
  cutoff <- as.Date(first - 28)
  dt[, date := substr(date, 1, 10)]
  dt[, date := as.Date(date, "%d/%m/%Y")]
  dt <- dt[date >= cutoff]
  
  #merge in municipality metadata
  dt <- merge(dt, muni.dt, by = "id_station")
  
  #Take the average temperature by day, for each station
  dt <- dt[, lapply(.SD, mean, na.rm = T), .SDcols = "temp_mean_day", by = c("date", "adm2_id", "adm2_name")] #min, max, mean
  #dt <- dt[, lapply(.SD, function(x){max(x, na.rm=T) - min(x, na.rm=T)}), .SDcols = "temp_mean_day", by = c("date", "adm2_id", "adm2_name")]
  
  #generating lagged mean temperature
  adm2s <- unique(dt$adm2_id)
  for (n in c(2,3,14,28)) {
    lag_coords_3 <- function(i, data){
      d <- data[adm2_id == i,]
      d[, paste0("temp_mean_",n,"day") := sapply(1:nrow(d), mean_lag_3, data = d, var = "temp_mean_day", lag = n)] ##
      return(d)
    }
  dt <- lapply(adm2s, lag_coords_3, data = dt) %>% rbindlist
  }
  #remove rows with previous year
  dt <- dt[date >= first]
  
  #generate full time series
  all_days <- seq.Date(from = as.Date(paste0(year, "-01-01")), to = as.Date(paste0(year, "-12-31")), by = "day")
  some_days <- dt[adm2_id== 9010, unique(date)]
  setdiff(some_days, all_days)
  
  #bind on to the end of main dt
  print(year)
  dt.all <- rbind(dt.all, dt)
}

#Clean and save
write.csv(dt.all, paste0(save.dir, "temp_mean_thesis.csv"), row.names = F)

##################################################################
#######pop prep
##################################################################
mex.city.munis <- c(9002, 9014, 9003, 9004, 9015, 9005, 9006, 9007, 9008, 9016, 9009, 9011, 9012, 9017, 9013, 9010)
dt.all <- data.table()
for(year in 1995:2016){
  dt <- fread(paste0(pop.dir, "lights_", year, ".csv"))
  dt <- dt[adm2_id %in% mex.city.munis,]
  dt.all <- rbind(dt.all, dt, fill = T)
}
dt.all <- dt.all[,.(population, adm2_id, year_id)]
write.csv(dt.all, paste0(save.dir, "pop.csv"), row.names = F)

##################################################################
#######pollutants prep
##################################################################
dt.all <- data.table()
for(year in 1998:2016){
  #read in data and subset to pollutant parameter
  dt <- fread(paste0(main.dir, "pollutants/contaminantes_", year, ".csv"))
  setnames(dt, c(colnames(dt)[grep("station", colnames(dt))], colnames(dt)[grep("parameter", colnames(dt))], "value"),
           c("id_station", "parameter", "o3_mean_day")) #min, max, mean
  dt2 <- fread(paste0(main.dir, "pollutants/contaminantes_", paste(year - 1), ".csv"))
  setnames(dt2, c(colnames(dt2)[grep("station", colnames(dt2))], colnames(dt2)[grep("parameter", colnames(dt2))], "value"),
           c("id_station", "parameter", "o3_mean_day")) #min, max, mean
  dt <- rbind(dt, dt2)
  rm(dt2)
  dt <- dt[parameter == "O3"] #PM10, PM2.5, O3
  
  #Read in year previous
  first <- as.Date(paste0("01/", as.numeric(year)), format = "%j/%Y")
  cutoff <- as.Date(first - 28)
  dt[, date := substr(date, 1, 10)]
  dt[, date := as.Date(date, "%d/%m/%Y")]
  dt <- dt[date >= cutoff]
  
  #merge in municipality metadata
  dt <- merge(dt, muni.dt, by = "id_station")
  
  #Take the average pollutant by day, for each station
  dt <- dt[, lapply(.SD, mean, na.rm =T), .SDcols = "o3_mean_day", by = c("date", "adm2_id", "adm2_name")]
  
  #generating lagged day mean pollutant concentration
  adm2s <- unique(dt$adm2_id)
  for (n in c(2,3,14,28)) {
    lag_coords_3 <- function(i, data){
      d <- data[adm2_id == i,]
      d[, paste0("o3_mean_",n,"day") := sapply(1:nrow(d), mean_lag_3, data = d, var = "o3_mean_day", lag = n)] ##
      return(d)
    }
    dt <- lapply(adm2s, lag_coords_3, data = dt) %>% rbindlist
  }
  #remove rows with previous year
  dt <- dt[date >= first]
  
  #generate full time series
  all_days <- seq.Date(from = as.Date(paste0(year, "-01-01")), to = as.Date(paste0(year, "-12-31")), by = "day")
  some_days <- dt[adm2_id== 9014, unique(date)]
  setdiff(some_days, all_days)
  
  #bind on to the end of main dt
  dt.all <- rbind(dt.all, dt)
  print(year)
}

#Clean and save
write.csv(dt.all, paste0(save.dir, "ozone_mean_thesis.csv"), row.names = F)

################################################################
#Merge them all together and save
predictors <- c("temp_min","temp_max","temp_range", "ozone", "pm10", "pm25", "so2", "humidity", "wind_speed") # will need to add population into this
predictors <- c("ozone", "pm10", "pm25")
#start with temperature
dt <- fread(paste0(save.dir, "temp_mean_thesis.csv"))
dt <- dt[, date := as.Date(date)]
for(p in predictors){
  dt.temp <- fread(paste0(save.dir, p, "_mean_thesis.csv"))
  dt.temp[, date := as.Date(date)]
  dt <- merge(dt, dt.temp, by = c("date", "adm2_id", "adm2_name"), all = T)
}
dt <- dt[!is.na(date)]
dt[, year_id := as.numeric(substr(date, 1, 4))]
dt.pop <- fread(paste0(save.dir, "pop.csv"))
dt <- merge(dt, dt.pop, by = c("adm2_id", "year_id"), all.x = T)
write.csv(dt, paste0(save.dir, "mex_preds_thesis.csv"), row.names = F)


##############################################################
#####cod data prep
##############################################################
muni.map <- fread(paste0(meta.dir, "cod_muni_map.csv"))
setnames(muni.map, c("CVE_ENT", "CVE_MUN"), c("ENT_OCURR", "MUN_OCURR"))
muni.map <- muni.map[!is.na(adm2_id)]
icd.map <- fread(paste0(meta.dir, "icd_map.csv"))

#broad icd codes
resp.codes <- paste0("J", 00:99)
cardio.codes <- paste0("I", 00:99)
ihd.codes <- paste0("I", 20:25)
stroke.codes <- paste0("I", 60:69)

dt.all <- data.table()
#Begin loop
for(year in 1998:2016){
  #Read in appropriate cod dataset
  ye <- substr(year, 3, 4)
  if(year < 2010){
    dt <- read.dbf(paste0(main.dir, "cod/defunciones_base_datos_", year, "_dbf/DEFUN", ye, ".dbf")) %>% as.data.table
  } else { dt <- read.dbf(paste0(main.dir, "cod/defunciones_2010_2016/DEFUN", ye, ".dbf")) %>% as.data.table }
  dt <- merge(dt, muni.map, by = c("ENT_OCURR", "MUN_OCURR"))
  
  #Generate date and age and sex variables
  dt[, date := as.Date(paste(DIA_OCURR, MES_OCURR, ANIO_OCUR, sep = "/"), "%e/%m/%Y")]
  dt[, birth_date := as.Date(paste(DIA_NACIM, MES_NACIM, ANIO_NACIM, sep = "/"), "%e/%m/%Y")]
  dt[, age := as.numeric(date - birth_date)/365]
  dt[, age := round(age)]
  dt[, age_bin := cut(age, seq(0,120,15), c(1:8))]
  dt[, sex_id := ifelse(SEXO == 1, 1, 2)]
  setnames(dt, c("ESCOLARIDA"), c("education"))
  
  #Merge on icd mapping
  dt[, all_cause_deaths := 1]
  dt[grepl(paste(cardio.codes, collapse = "|"), dt$CAUSA_DEF), cardio_deaths := 1]
  dt[is.na(cardio_deaths), cardio_deaths := 0]
  dt[grepl(paste(resp.codes, collapse = "|"), dt$CAUSA_DEF), resp_deaths := 1]
  dt[is.na(resp_deaths), resp_deaths := 0]
  dt[grepl(paste(ihd.codes, collapse = "|"), dt$CAUSA_DEF), ihd_deaths := 1]
  dt[is.na(ihd_deaths), ihd_deaths := 0]
  dt[grepl(paste(stroke.codes, collapse = "|"), dt$CAUSA_DEF), stroke_deaths := 1]
  dt[is.na(stroke_deaths), stroke_deaths := 0]
  
  #collapse to amd2, day level
  dt <- dt[, lapply(.SD, sum), .SDcols = c("all_cause_deaths", "cardio_deaths", "resp_deaths", "ihd_deaths", "stroke_deaths"), by = c("adm2_id", "date")]
  #dt <- dt[, lapply(.SD, sum), .SDcols = c("all_cause_deaths", "cardio_deaths", "resp_deaths", "ihd_deaths", "stroke_deaths"), by = c("adm2_id", "date", "sex_id", "age_bin")]
  dt.all <- rbind(dt.all, dt)
  print(year)
}

###Read in predictor variables, format, and merge####################
dt.preds <- fread(paste0(save.dir, "mex_preds_thesis.csv"))
dt.preds[, date := as.Date(date)]
dt <- merge(dt.all, dt.preds, by= c("adm2_id", "date"), all.y = T)
dt <- dt[!is.na(date)]
dt <- dt[date > as.Date("1997-12-31") & date < as.Date("2017-01-01")]
dt[is.na(stroke_deaths), stroke_deaths := 0]
dt[is.na(ihd_deaths), ihd_deaths := 0]
dt[is.na(resp_deaths), resp_deaths := 0]
dt[is.na(cardio_deaths), cardio_deaths := 0]
dt[is.na(all_cause_deaths), all_cause_deaths := 0]

#Create season variable
get_season <- function(dates) {
  #set dates that delineate seasons
  WS <- as.Date("2012-12-15", format = "%Y-%m-%d") # Winter Solstice
  SE <- as.Date("2012-3-15",  format = "%Y-%m-%d") # Spring Equinox
  SS <- as.Date("2012-6-15",  format = "%Y-%m-%d") # Summer Solstice
  FE <- as.Date("2012-9-15",  format = "%Y-%m-%d") # Fall Equinox
  
  # Convert dates from any year to 2012 dates bc it's a leap year
  d <- as.Date(strftime(dates, format="2012-%m-%d"))
  
  ifelse (d >= WS | d < SE, "Winter",
          ifelse (d >= SE & d < SS, "Spring",
                  ifelse (d >= SS & d < FE, "Summer", "Fall")))
}
dt[, season := get_season(date)]

#Add trend variable
dt <- dt[order(adm2_id, date)]
dt[, trend := sequence(.N), by = adm2_id]
dt[, day_o_week := weekdays(date)]
write.csv(dt, paste0(model.dir, "mex_model_thesis.csv"), row.names = F)


###SCRAP########################################################
# dt <- fread(paste0(main.dir, "cod/all_causes_mexcity.csv"))
# dt <- dt[!is.na(mmt)]
# dt <- dt[,.(adm2_id_res, adm2_name_res, n_cvd, n_resp, n__ri, n_all_cause, date)]
# dt[, date := as.Date(date, "%d%b%Y")]
# setnames(dt, "adm2_id_res", "adm2_id")
# 
# #Read in predictor variables, format, and merge
# dt.preds <- fread(paste0(save.dir, "all_preds.csv"))
# dt.preds[, date := as.Date(date)]
# dt <- merge(dt, dt.preds, by= c("adm2_id", "date"), all.x = T)
# dt[is.na(n_cvd), n_cvd := 0]
# dt[is.na(n_resp), n_resp := 0]
# dt[is.na(n__ri), n__ri := 0]
# write.csv(dt, paste0(model.dir, "/mex_city_full.csv"), row.names = F)
