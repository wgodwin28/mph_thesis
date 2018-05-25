#Clean and prep meteorology, pollution variables
rm(list=ls())
#set up smart sensing of os
if(Sys.info()[1]=='Windows'){
  j = "J:/"
} else{
  j = '/home/j/'
}

# load packages, install if missing 
pacman::p_load(data.table, fst, ggplot2, parallel, magrittr, lavaan, lavaanPlot)

#set settings
main.dir <- paste0(j, "temp/wgodwin/thesis/data/01_raw/")
save.dir <- paste0(j, "temp/wgodwin/thesis/data/02_prepped/")
model.dir <- paste0(j, "temp/wgodwin/thesis/data/03_analysis/")
meta.dir <- paste0(j, "temp/wgodwin/thesis/metadata/")
pop.dir <- paste0(j, "WORK/05_risk/risks/temperature/data/exp/pop/mex/all_age/")
muni.dt <- fread(paste0(meta.dir, "station_muni_map.csv"))
day <- F
######################################Function to calculate moving average over 30 day lag############################
mean_lag_3 <- function(i, data, var) {
  is.near <- as.numeric(data$date[i] - data$date) >= 0 & as.numeric(data$date[i] - data$date) < 3
  #mean(data$paste(var)[is.near], na.rm = T)
  data[, max(get(var)[is.near], na.rm = T)]
}

##################################################################
#######temperature prep
##################################################################
dt.all <- data.table()
for(year in 1995:2015){
  #read in data and subset to temperature parameter
  dt <- fread(paste0(main.dir, "meteorology/meteorología_", year, ".csv"))
  setnames(dt, c(colnames(dt)[grep("station", colnames(dt))], colnames(dt)[grep("parameter", colnames(dt))], "value"),
           c("id_station", "parameter", "temp_max_day")) #min, max, mean
  dt <- dt[parameter == "TMP"] #RH, TMP, WSP
  
  #reformat date
  dt[, date := substr(date, 1, 10)]
  dt[, date := as.Date(date, "%d/%m/%Y")]
  
  #merge in municipality metadata
  dt <- merge(dt, muni.dt, by = "id_station")
  
  #Take the average temperature by day, for each station
  #if(day){
    dt <- dt[, lapply(.SD, max, na.rm = T), .SDcols = "temp_max_day", by = c("date", "adm2_id", "adm2_name")] #min, max, mean
    #dt <- dt[, lapply(.SD, function(x){max(x, na.rm=T) - min(x, na.rm=T)}), .SDcols = "temp_range_day", by = c("date", "adm2_id", "adm2_name")]
  #}else{
    #generating previous 30 day mean temperature
    adm2s <- unique(dt$adm2_id)
    lag_coords_3 <- function(i, data){
      d <- data[adm2_id == i,]
      d[, temp_range_3day := sapply(1:nrow(d), mean_lag_3, data = d, var = "temp_max_day")] ##
      return(d)
    }
    dt <- lapply(adm2s, lag_coords_3, data = dt) %>% rbindlist
  #}
  #bind on to the end of main dt
  print(year)
  dt.all <- rbind(dt.all, dt)
}

#Clean and save
write.csv(dt.all, paste0(save.dir, "temp_max.csv"), row.names = F)

##################################################################
#######pop prep
##################################################################
mex.city.munis <- c(9002, 9014, 9003, 9004, 9015, 9005, 9006, 9007, 9008, 9016, 9009, 9011, 9012, 9017, 9013, 9010)
dt.all <- data.table()
for(year in 1995:2015){
  dt <- fread(paste0(pop.dir, "lights_", year, ".csv"))
  dt <- dt[adm2_id %in% mex.city.munis,]
  dt.all <- rbind(dt.all, dt, fill = T)
}
dt <- dt[,.(population, adm2_id, year_id)]
write.csv(dt, paste0(save.dir, "pop.csv"), row.names = F)

##################################################################
#######pollutants prep
##################################################################
dt.all <- data.table()
for(year in 1995:2015){
  #read in data and subset to temperature parameter
  dt <- fread(paste0(main.dir, "pollutants/contaminantes_", year, ".csv"))
  setnames(dt, c(colnames(dt)[grep("station", colnames(dt))], colnames(dt)[grep("parameter", colnames(dt))], "value"), c("id_station", "parameter", "pm25"))
  dt <- dt[parameter == "PM2.5"]
  
  #reformat date
  dt[, date := substr(date, 1, 10)]
  dt[, date := as.Date(date, "%d/%m/%Y")]
  
  #merge in municipality metadata
  dt <- merge(dt, muni.dt, by = "id_station")
  
  #Take the average temperature by day, for each station
  dt <- dt[, lapply(.SD, mean, na.rm =T), .SDcols = "pm25", by = c("date", "adm2_id", "adm2_name")]
  
  #bind on to the end of main dt
  dt.all <- rbind(dt.all, dt)
}

#Clean and save
write.csv(dt.all, paste0(save.dir, "pm25.csv"), row.names = F)

################################################################
#Merge them all together and save
predictors <- c("temp_min","temp_max","temp_range", "ozone", "pm10", "pm25", "so2", "humidity", "wind_speed") # will need to add population into this

#start with temperature
dt <- fread(paste0(save.dir, "temp_mean.csv"))
dt <- dt[, date := as.Date(date)]
for(p in predictors){
  dt.temp <- fread(paste0(save.dir, p, ".csv"))
  dt.temp[, date := as.Date(date)]
  dt <- merge(dt, dt.temp, by = c("date", "adm2_id", "adm2_name"), all.x = T)
}
dt <- dt[!is.na(date)]
dt[, year_id := as.numeric(substr(date, 1, 4))]
dt.pop <- fread(paste0(save.dir, "pop.csv"))
dt <- merge(dt, dt.pop, by = c("adm2_id", "year_id"), all.x = T)
write.csv(dt, paste0(model.dir, "mex_city_full.csv"), row.names = F)

##############################################################
#####run model
##############################################################
dt <- fread(paste0(j, "temp/wgodwin/thesis/data/03_analysis/mex_city_full.csv"))
model <- 'pm25 ~ temperature + humidity + wind_speed'

#model 1
sem_mod_1 <- "temperature  =~ temp_mean_day + temp_min_day + temp_max_day
              ozone ~ temperature + humidity + wind_speed"

sem_mod_1 <- "ozone ~ temp_min_day + humidity + wind_speed"

fit <- sem(model = sem_mod_1, data = dt)
summary(fit)
lavaanPlot(model = fit)
dt[temp_min_3day< -100000, temp_min_3day := NA]
dt[temp_min_day< -100000, temp_min_day := NA]
dt[temp_max_day< -100000, temp_max_day := NA]
dt[temp_max_3day< -100000, temp_max_3day := NA]
dt[temp_max_3day> 100000, temp_max_3day := NA]
dt[temp_max_day> 100000, temp_max_day := NA]
dt[temp_min_day> 100000, temp_min_day := NA]
dt[temp_min_3day> 100000, temp_min_3day := NA]
dt[temp_max_3day> 100000, temp_max_3day := NA]

