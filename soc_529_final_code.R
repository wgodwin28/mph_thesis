#Clean and prep meteorology, pollution variables
rm(list=ls())
#set up smart sensing of os
if(Sys.info()[1]=='Windows'){
  j = "J:/"
} else{
  j = '/home/j/'
}

# load packages, install if missing 
pacman::p_load(data.table, fst, ggplot2, parallel, magrittr, lavaan, lavaanPlot, semPlot, xtable, Hmisc)

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
dt <- fread("C:/Users/wgodwin/Desktop/csss_526/exercises/mex_city_full.csv")
dt2 <- dt[,.(temp_mean_day, temp_mean_3day, temp_min_day, temp_min_3day, temp_max_day, temp_max_3day, humidity, wind_speed, ozone, pm25)]
mcor <- round(cor(dt2, use = "complete.obs"), 3)
symnum(corre.matr)
lower.tri(dt2, diag = F)
upper<-mcor
upper[upper.tri(mcor)]<-""
upper<-as.data.frame(upper)
pdf("C:/Users/wgodwin/Desktop/csss_526/exercises/mcor.pdf")
print(xtable(upper), type="html")
dev.off()
corstars(dt2, result="html")
corrplot(mcor, method = "number")

dt2 <- melt(dt2, measure.vars = colnames(dt2))
pp <- ggplot(aes(x = value), data = dt2) + facet_wrap(~variable, scales="free") + 
      geom_histogram() + stat_function(fun = function(x, mean, sd, n){ 
        n * dnorm(x = value, mean = mean, sd = sd)}, 
        args = list(mean = mean(dt2$value), sd = sd(dt2$value), n = 100))
#####model 1######
sem_mod_1 <- "temperature  =~ temp_mean_day + temp_min_day + temp_max_day
              ozone ~ temperature + humidity + wind_speed"

fit1 <- cfa(model = sem_mod_1, data = dt)
summary(fit1, fit.measures = T)
lavaanPlot(model = fit1,  stand = F, coefs = T, sig = .005, labels = c("ozone", "Ozone", "Temperature Maximum"))
semPaths(fit1)

#####model 2######
sem_mod_2 <- "temperature  =~ temp_mean_3day + temp_min_3day + temp_max_3day
              ozone ~ temperature + humidity + wind_speed"

fit2 <- cfa(model = sem_mod_2, data = dt)
summary(fit2, fit.measures = T)
lavaanPlot(model = fit2,  coefs = T)


#####model 3######
sem_mod_3 <- 
"temperature_daily  =~ temp_mean_day + temp_min_day + temp_max_day
temperature_3day    =~ temp_mean_3day + temp_min_3day + temp_max_3day
ozone ~ temperature_daily + temperature_3day + humidity + wind_speed"

fit3 <- cfa(model = sem_mod_3, data = dt)
summary(fit3, fit.measures = T)
lavaanPlot(model = fit3,  coefs = T)

#####model 4######
sem_mod_4 <- "temperature  =~ temp_mean_day + temp_min_day + temp_max_day
              pm25 ~ temperature + humidity + wind_speed"

fit4 <- cfa(model = sem_mod_4, data = dt)
summary(fit4, fit.measures = T)
lavaanPlot(model = fit4,  coefs = T)


#####model 5######
sem_mod_5 <- "temperature  =~ temp_mean_3day + temp_min_3day + temp_max_3day
              pm25 ~ temperature + humidity + wind_speed"

fit5 <- cfa(model = sem_mod_5, data = dt)
summary(fit5, fit.measures = T)
lavaanPlot(model = fit5,  coefs = T)


#####model 6######
sem_mod_6 <- 
  "temperature_daily  =~ temp_mean_day + temp_min_day + temp_max_day
  temperature_3day    =~ temp_mean_3day + temp_min_3day + temp_max_3day
  pm25 ~ temperature_daily + temperature_3day + humidity + wind_speed"

fit6 <- cfa(model = sem_mod_6, data = dt)
summary(fit6, fit.measures = T)
lavaanPlot(model = fit6,  coefs = T, stand = F)


dt[temp_min_3day< -100000, temp_min_3day := NA]
dt[temp_min_day< -100000, temp_min_day := NA]
dt[temp_max_day< -100000, temp_max_day := NA]
dt[temp_max_3day< -100000, temp_max_3day := NA]
dt[temp_max_3day> 100000, temp_max_3day := NA]
dt[temp_max_day> 100000, temp_max_day := NA]
dt[temp_min_day> 100000, temp_min_day := NA]
dt[temp_min_3day> 100000, temp_min_3day := NA]
dt[temp_max_3day> 100000, temp_max_3day := NA]

##################################################################################
#######################THESIS DATA################################################
##################################################################################
library(dlnm) ; library(mvmeta) ; library(splines) ; library(tsModel) ; library(data.table) ; library(zoo) ; library(mgcv) ; library(itsadug)

#Directories
in.dir <- paste0(j, "temp/wgodwin/thesis/data/03_analysis/")

###Some initial exploratory plots of the exposure and response variables
#Read in data
dt <- fread("C:/Users/wgodwin/Desktop/csss_526/exercises/mex_city_full.csv")
dt[, date := as.Date(date)]
dt[, month := months(date, abbreviate=T)]

corstars <-function(x, method=c("pearson", "spearman"), removeTriangle=c("upper", "lower"),
                    result=c("none", "html", "latex")){
  #Compute correlation matrix
  require(Hmisc)
  x <- as.matrix(x)
  correlation_matrix<-rcorr(x, type=method[1])
  R <- correlation_matrix$r # Matrix of correlation coeficients
  p <- correlation_matrix$P # Matrix of p-value 
  
  ## Define notions for significance levels; spacing is important.
  mystars <- ifelse(p < .0001, "****", ifelse(p < .001, "*** ", ifelse(p < .01, "**  ", ifelse(p < .05, "*   ", "    "))))
  
  ## trunctuate the correlation matrix to two decimal
  R <- format(round(cbind(rep(-1.11, ncol(x)), R), 2))[,-1]
  
  ## build a new matrix that includes the correlations with their apropriate stars
  Rnew <- matrix(paste(R, mystars, sep=""), ncol=ncol(x))
  diag(Rnew) <- paste(diag(R), " ", sep="")
  rownames(Rnew) <- colnames(x)
  colnames(Rnew) <- paste(colnames(x), "", sep="")
  
  ## remove upper triangle of correlation matrix
  if(removeTriangle[1]=="upper"){
    Rnew <- as.matrix(Rnew)
    Rnew[upper.tri(Rnew, diag = TRUE)] <- ""
    Rnew <- as.data.frame(Rnew)
  }
  
  ## remove lower triangle of correlation matrix
  else if(removeTriangle[1]=="lower"){
    Rnew <- as.matrix(Rnew)
    Rnew[lower.tri(Rnew, diag = TRUE)] <- ""
    Rnew <- as.data.frame(Rnew)
  }
  
  ## remove last column and return the correlation matrix
  Rnew <- cbind(Rnew[1:length(Rnew)-1])
  if (result[1]=="none") return(Rnew)
  else{
    if(result[1]=="html") print(xtable(Rnew), type="html")
    else print(xtable(Rnew), type="latex") 
  }
} 
