#Model
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
library(dlnm) ; library(mvmeta) ; library(splines) ; library(tsModel) ; library(data.table) ; library(zoo) ; library(mgcv) ; library(itsadug)

#Directories
in.dir <- paste0(j, "temp/wgodwin/thesis/data/03_analysis/")

###Some initial exploratory plots of the exposure and response variables
#Read in data
dt <- fread(paste0(in.dir, "/mex_model_thesis.csv"))
dt[, date := as.Date(date)]
dt[, month := months(date, abbreviate=T)]

#Temperature and pollutant collapse
plot.dt <- dt[, lapply(.SD, mean, na.rm = T), by = "month", .SDcols = c("temp_mean_day", "o3_mean_day", "pm25_mean_day")]
plot.dt <- melt(plot.dt, id.vars = "month")
plot.dt[]

#plot
ggplot(plot.dt, aes(month, value, group = variable, color = variable)) +
  geom_point() + geom_line() + scale_x_discrete(limits = month.abb)

#Mortality collapse
plot.dt <- dt[, lapply(.SD, sum, na.rm = T), by = c("adm2_id", "month"), .SDcols = c("population", "all_cause_deaths")]
temp <- dt[, lapply(.SD, sum, na.rm = T), by = c("month"), .SDcols = c("population", "all_cause_deaths")]
temp[, adm2_id := 1111]
plot.dt <- rbind(plot.dt, temp)
plot.dt[, death_rate := all_cause_deaths/population * 365.25] #converts denominator of person-days to person-years
plot.dt[, death_rate := death_rate * 1000] #converts the rate to per 1000 person-years

#plot
ggplot(plot.dt, aes(month, death_rate, group = as.factor(adm2_id), color = as.factor(adm2_id))) +
  geom_point() + geom_line() + scale_x_discrete(limits = month.abb)

#Summary statistics for a table similar to K's paper
dt[,sum(all_cause_deaths)]
summary(dt$temp_mean_day)
quantile(dt$temp_mean_day, c(.05, .95), na.rm = T)
summary(dt$o3_mean_day)
quantile(dt$o3_mean_day, c(.05, .95), na.rm = T)
summary(dt$pm10_mean_day)
quantile(dt$pm10_mean_day, c(.05, .95), na.rm = T)
summary(dt$pm25_mean_day)
quantile(dt$pm25_mean_day, c(.05, .95), na.rm = T)


##############GAMs##########################
dt[, year_id := as.numeric(substr(date, 1, 4))]

#start with something basic (just temperature and pollution variable)
gamO3.bas <- gam(all_cause_deaths ~ s(temp_mean_2day, bs= "cs", k = 10) + s(o3_mean_2day, bs= "cs", k = 10) + offset(log(population)),
             family=poisson, data=dt)
plot(gamO3.bas, select = 2, shade = T, rug = T) # residuals = T

#move to adding trend variables- seasonal variables significant and in expected direction (summer and spring with lowest mortality)
  #heat effects look to be strongest. Why?
gamO3.bas2 <- gam(all_cause_deaths ~ s(temp_mean_2day, bs= "cs", k = 10) + s(o3_mean_2day, bs= "cs", k = 10) + as.factor(season) +
                  s(trend,bs="cs", k=67, fx=T) + offset(log(population)),
                 family=poisson, data=dt)
plot(gamO3.bas2, select = 2, shade = T, rug = T) # residuals = T

#testing out random effects on municipality
dt[, adm2_id := as.factor(adm2_id)]
gamO3.re <- gam(all_cause_deaths ~ s(temp_mean_2day, bs= "cs", k = 10) + s(o3_mean_2day, bs= "cs", k = 10) + as.factor(season) +
                    s(trend, bs="cs", k=67, fx=T) + s(adm2_id, bs = "re") + offset(log(population)),
                  family=poisson, data=dt)
plot(gamO3.re, select = 1)

  #look at predictions and random effects
  plot_smooth(gamO3.re, view = "temp_mean_2day", rm.ranef = T) # w/o random effects
  plot_smooth(gamO3.re, view = "temp_mean_2day", cond = list(adm2_id = 9007), rug = F, col = 'red')
  plot_smooth(gamO3.re, view = "temp_mean_2day", cond = list(adm2_id = 9010), add = T, col = 'blue', xpd = T)
  plot_smooth(gamO3.re, view = "temp_mean_2day", plot_all = "adm2_id")


#Make ozone and temperature a tensor product
gamO3.ten <- gam(all_cause_deaths ~  te(temp_mean_2day, o3_mean_2day, k=c(4,5), fx = T) + as.factor(season) +
                   s(trend,bs="cs",k=67,fx=T) + s(adm2_id, bs = "re") + offset(log(population)),
                 family=poisson, data=dt)
plot(gamO3.ten, select = 1, pers = T,ticktype="detailed", theta = -35, too.far = .15, se = T, col = "orange2",
     main = "Log relative risk of mortality", xlab = "2-days average temperature [°C]", ylab = "2-day average ozone [µg/m³]") # residuals = T
plot(gamO3.ten, all.terms = T)

#Make ozone and temperature a tensor product-decrease degrees of freedom on temp and ozone
gamO3.ten.sm <- gam(all_cause_deaths ~  te(temp_mean_2day, o3_mean_2day, k=c(3,4), fx = T) + as.factor(season) +
                      s(trend,bs="cs",k=67,fx=T) + s(adm2_id, bs = "re") + offset(log(population)),
                    family=poisson, data=dt)
plot(gamO3.ten.sm, select = 1, pers = T, theta = -35, se=T, ticktype="detailed") # residuals = T

#Make ozone and temperature a tensor product-try month trend instead of season**************************BEST MODEL************************
gamO3.ten.mon <- gam(all_cause_deaths ~  te(temp_mean_2day, o3_mean_2day, k=c(4,5), fx = T) + as.factor(month) +
                   s(trend,bs="cs",k=67,fx=T) + s(adm2_id, bs = "re") + offset(log(population)),
                 family=poisson, data=dt)
plot(gamO3.ten.mon, select = 1, pers = T,ticktype="detailed", theta = -35, too.far = .15, se = T, col = "orange2", 
       main = "Log relative risk of mortality", xlab = "2-days average temperature [°C]", ylab = "2-day average ozone [µg/m³]") # residuals = T
plot(gamO3.ten.mon, all.terms = T)

#Make ozone and temperature a tensor product-include month and season
gamO3.ten.monsea <- gam(all_cause_deaths ~  te(temp_mean_2day, o3_mean_2day, k=c(4,5), fx = T) + as.factor(month) + as.factor(season) +
                       s(trend,bs="cs",k=67,fx=T) + s(adm2_id, bs = "re") + offset(log(population)),
                     family=poisson, data=dt)
plot(gamO3.ten.monsea, select = 1, pers = T,ticktype="detailed", theta = -35, too.far = .15, se = T, col = "orange2", 
     main = "Log relative risk of mortality", xlab = "2-days average temperature [°C]", ylab = "2-day average ozone [µg/m³]") # residuals = T
plot(gamO3.ten.monsea, all.terms = T)

#Make ozone and temperature a tensor product-only include month and reduce knots on trend
gamO3.ten.monsm <- gam(all_cause_deaths ~  te(temp_mean_2day, o3_mean_2day, k=c(4,5), fx = T) + as.factor(month) +
                          s(trend,bs="cs",k=40,fx=T) + s(adm2_id, bs = "re") + offset(log(population)),
                        family=poisson, data=dt)
plot(gamO3.ten.monsm, select = 1, pers = T,ticktype="detailed", theta = -35, too.far = .15, se = T, col = "orange2", 
     main = "Log relative risk of mortality", xlab = "2-days average temperature [°C]", ylab = "2-day average ozone [µg/m³]") # residuals = T
plot(gamO3.ten.monsm, all.terms = T)

#Make ozone and temperature a tensor product-only include month and increase knots on trend
gamO3.ten.monlesm <- gam(all_cause_deaths ~  te(temp_mean_2day, o3_mean_2day, k=c(4,5), fx = T) + as.factor(month) +
                         s(trend,bs="cs",k=80,fx=T) + s(adm2_id, bs = "re") + offset(log(population)),
                       family=poisson, data=dt)
plot(gamO3.ten.monlesm, select = 1, pers = T,ticktype="detailed", theta = -35, too.far = .15, se = T, col = "orange2", 
     main = "Log relative risk of mortality", xlab = "2-days average temperature [°C]", ylab = "2-day average ozone [µg/m³]") # residuals = T
plot(gamO3.ten.monlesm, all.terms = T)

# #Switch ozone and temperature to be modeled without smoothing-compare with model that puts splines on them-this doesn't make sense to do since temp function isn't linear
# gamO3.nospl <- gam(all_cause_deaths ~ temp_mean_2day + o3_mean_2day + as.factor(season) +
#                     s(trend,bs="cs",k=67,fx=T) + offset(log(population)),
#                   family=poisson, data=dt)
# plot(gamO3.nospl, shade = T, rug = T) # residuals = T

gamO3 <- gam(all_cause_deaths ~ te(trend, k=67, fx=T) + as.factor(season) + as.factor(day_o_week) + as.factor(season)
                 te(temp_mean_2day, o3_mean_2day, k=c(4,5), fx = T) + offset(log(population)), 
                 family=poisson, data=dt)

gamO3 <- gam(all_cause_deaths ~ te(trend, k=40, fx=T) + as.factor(season) + as.factor(day_o_week) + 
               te(temp_mean_2day, o3_mean_2day, k=c(4,5), fx = T) + offset(log(population)), 
             family=poisson, data=dt)
gamO3.year <- gam(all_cause_deaths ~ as.factor(year_id) + as.factor(season) + as.factor(day_o_week) + 
               te(temp_mean_2day, o3_mean_2day, k=c(4,5), fx = T) + offset(log(population)), 
             family=poisson, data=dt)
#add population. Vary df for trend variable and include season for sensitivity analysis
#Seperating long term from short term. Patrick Kinney
#play around with number of cats for interaction
#Adjust dof in the te product and play around with different values
#use quasipossion instead poisson family
plot(gamO3, select = 2, scheme = 2, theta = -35, ticktype="detailed", pers = T)
plot(gamO3.year, select = 2, scheme = 2, theta = -35, ticktype="detailed", pers = T)

#plot 3 types (3-D, trend plot, contour)
plot(gamO3, scheme = 2, theta = -35, ticktype="detailed", pers = F)
par(mfrow=c(1,2))
axis(4)
plot(gamO3, select = 1, scheme = 2, theta = -35, ticktype="detailed")
vis.gam(gamO3, view = c("temp_mean_2day", "o3_mean_2day"), type = "response", se = 1)
