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
library(dlnm) ; library(mvmeta) ; library(splines) ; library(tsModel) ; library(zoo) ; library(mgcv) ; library(itsadug) ; library(hexbin) ; library(magrittr) ; library(ggplot2) ; library(kableExtra) ; library(knitr) ; library(RColorBrewer) ; library(data.table)

#Directories
in.dir <- paste0(j, "temp/wgodwin/thesis/data/03_analysis/")
pred.dir <- paste0(j, "temp/wgodwin/thesis/data/02_prepped/")

###Some initial exploratory plots of the exposure and response variables
#Read in data
dt <- fread(paste0(in.dir, "/mex_model_thesis.csv"))
#dt <- fread(paste0(in.dir, "/mex_model_thesis_age65.csv"))
dt[, adm2_id := as.factor(adm2_id)]
dt[, date := as.Date(date)]
dt[, month := months(date, abbreviate=T)]
dt[, year_id := as.numeric(substr(date, 1, 4))]
days <- seq(as.Date("1998-01-01"), as.Date("2016-12-31"), by="days")
t <- data.table(date = days, trend2 = 1:length(days))
dt <- merge(dt, t, by = "date", all.x = T)
dt[, adm2_id := as.integer(as.character(adm2_id))]

#Merge on range, min, max
predictors <- c("range", "max", "min")
for(p in predictors){
  t <- fread(paste0(pred.dir, "temp_", p, ".csv"))
  t$adm2_name <- NULL
  t[, date := as.Date(date)]
  dt <- merge(dt, t, by = c("date", "adm2_id"), all.x = T)
}
dt[temp_max_day==-Inf, temp_max_day := NA]
dt[temp_min_day==-Inf, temp_min_day := NA]
dt[temp_range_day==-Inf, temp_range_day := NA]
dt[, adm2_id := as.factor(adm2_id)]

#Scatters
ls <- loess(temp_mean_day ~ o3_mean_day, data = dt)
pr.loess <- predict(ls)

ggplot(data = dt[year_id == 2016], aes(temp_mean_day, o3_mean_day)) + geom_point() + geom_smooth(method = "lm")
ggplot(data = dt[year_id == 2003], aes(temp_mean_day, pm25_mean_day)) + geom_point() + geom_smooth(method = "loess")
ggplot(data = dt, aes(temp_mean_day, temp_mean_2day)) + geom_point() + geom_smooth(method = "lm")
ggplot(data = dt, aes(o3_mean_day, o3_mean_2day)) + geom_point() + geom_smooth(method = "lm")

dt[, deaths_prop := (all_cause_deaths/population)*1000]
ggplot(data = dt, aes(o3_mean_day, deaths_prop)) + geom_point() + geom_smooth(method = "lm")
plot(dt$temp_mean_day, dt$o3_mean_day)
!is.na(dt$temp_mean_day)
lines(pr.loess~!is.na(dt$temp_mean_day), col="red", lwd=2)

rf <- colorRampPalette(rev(brewer.pal(11,'Spectral')))
hexbinplot(o3_mean_day ~ temp_mean_day, data=dt, colramp=rf, ylab="Mean Daily Ozone Concentration", xlab = "Mean Daily Temperature")
hexbinplot(pm25_mean_day ~ temp_mean_day, data=dt[pm25_mean_day < 130], colramp=rf, ylab="Mean Daily PM 2.5 Concentration", xlab = "Mean Daily Temperature")
hexbinplot(cardio_deaths ~ pm25_mean_2day, data=dt, colramp=rf)

#ME linear regression
mod <- lm(all_cause_deaths ~ o3_mean_day + temp_mean_day, data = dt)
mod.me <- lmer(all_cause_deaths ~ o3_mean_day + temp_mean_day + (1|adm2_id), data = dt)
mod.pos <- glmer(all_cause_deaths ~ o3_mean_day + (1|adm2_id), offset = log(population), data = dt, family = poisson)

#munis represented in muni-years
t <- dt[, lapply(.SD, sum), .SDcols = "all_cause_deaths", by=c("adm2_id","adm2_name", "year_id")]
t[, count := 1]
t[, lapply(.SD, sum), .SDcols = "count", by = c("adm2_id", "adm2_name")]

#Correlation matrix
cor.dt <- dt[,.(temp_mean_day, o3_mean_day, pm25_mean_day)]
setnames(cor.dt, c("temp_mean_day", "o3_mean_day", "pm25_mean_day"), c("Temperature", "Ozone", "PM 2.5"))
t <- cor(cor.dt, use = "complete.obs") %>% round(digits = 2)
t[lower.tri(t, diag=F)]<-""
#pdf("/home/j/temp/wgodwin/test.pdf")
kable(t, booktabs = T)
dev.off()
print(xtable(t))
cov(cor.dt, use = "complete.obs")
sjt.corr(cor.dt)

#Temperature and pollutant collapse
plot.dt <- dt[, lapply(.SD, mean, na.rm = T), by = c("season", "year_id"), .SDcols = c("temp_mean_day", "o3_mean_day", "pm25_mean_day")]
plot.dt <- melt(plot.dt, id.vars = c("season", "year_id"))
plot.dt[, day := 01]
plot.dt[season == "Winter", month := "Jan"]
plot.dt[season == "Spring", month := "Apr"]
plot.dt[season == "Summer", month := "Jul"]
plot.dt[season == "Fall", month := "Oct"]

plot.dt[, time := as.Date(paste0(day, "-", month, "-", year_id), format = "%d-%b-%Y")]
plot.dt[, time := as.Date(time, "%b-%Y")]
plot.dt[variable == "temp_mean_day", Variable := "Temperature (°C)"]
plot.dt[variable == "o3_mean_day", Variable := "Ozone (ppb)"]
plot.dt[variable == "pm25_mean_day", Variable := "PM 2.5 (ppb)"]
plot.dt[value > 55, value := NA]

#plot
ggplot(plot.dt, aes(time, value, group = Variable, color = Variable)) +
  geom_point() + geom_line() + theme_bw() + scale_x_date() +  xlab("Year") +
  theme(text = element_text(size=23)) + scale_color_manual(values=c("#56B4E9", "orange", "green"))
  #ylab("Degrees Celsius (temperature) or PPB (pollutants)")

#Mortality collapse
plot.dt <- copy(dt)
setnames(plot.dt, "all_cause_deaths", "deaths")
plot.dt <- plot.dt[, lapply(.SD, sum, na.rm = T), by = c("adm2_id", "month"), .SDcols = c("population", "deaths")]
temp <- plot.dt[, lapply(.SD, sum, na.rm = T), .SDcols = c("population", "deaths")]
temp[, adm2_id := 1111]
plot.dt <- rbind(plot.dt, temp)
plot.dt[, death_rate := (deaths/population) * 365.25] #converts denominator of person-days to person-years
plot.dt[, death_rate := death_rate * 1000] #converts the rate to per 1000 person-years

plot.dt <- copy(dt)
setnames(plot.dt, "all_cause_deaths", "deaths")
plot.dt <- plot.dt[, lapply(.SD, sum, na.rm = T), by = c("year_id", "adm2_id"), .SDcols = c("deaths")]
setkeyv(dt, c("adm2_id", "year_id"))
dt[ , i := .GRP, by = key(dt)]
pop <- dt[!duplicated(i), .(adm2_id, population, year_id)]
plot.dt <- merge(plot.dt, pop, by = c("adm2_id", "year_id"))
plot.dt[, lapply(.SD, sum, na.rm = T), .SDcols = c("deaths", "population")]

#plot
ggplot(plot.dt, aes(month, death_rate, group = as.factor(adm2_id), color = as.factor(adm2_id))) +
  geom_point() + geom_line() + scale_x_discrete(limits = month.abb) + theme_bw()

####Check on long term death trend###
plot.dt <- copy(dt)
#plot.dt <- plot.dt[adm2_id != 9015]
setnames(plot.dt, "all_cause_deaths", "deaths")
plot.dt <- plot.dt[, lapply(.SD, sum, na.rm = T), by = c("adm2_id", "season", "year_id"), .SDcols = c("population", "deaths")]
temp <- plot.dt[, lapply(.SD, sum, na.rm = T), by = c("season", "year_id"), .SDcols = c("population", "deaths")]
temp[, adm2_id := 1111]
plot.dt <- rbind(plot.dt, temp)
plot.dt[, death_rate := (deaths/population) * 365.25] #converts denominator of person-days to person-years
plot.dt[, death_rate := death_rate * 1000]
plot.dt[, day := 01]
plot.dt[season == "Winter", month := "Jan"]
plot.dt[season == "Spring", month := "Apr"]
plot.dt[season == "Summer", month := "Jul"]
plot.dt[season == "Fall", month := "Oct"]
plot.dt[, time := as.Date(paste0(day, "-", month, "-", year_id), format = "%d-%b-%Y")]
plot.dt[, time := as.Date(time, "%b-%Y")]

#plot
ggplot(plot.dt[adm2_id == 1111], aes(time, death_rate, group = as.factor(adm2_id), color = as.factor(adm2_id))) +
  geom_point() + geom_line() + scale_x_date()
ggplot(plot.dt[adm2_id == 1111], aes(time, death_rate, color = "red")) +
  geom_point() + geom_line() + scale_x_date() + xlab("Year") + ylab("Death Rate (per 1000)") +
  theme_bw() + theme(text = element_text(size=20))

#Summary statistics for a table similar to K's paper
dt[, sum(resp_deaths)]
dt[, lapply(.SD, sum), .SDcols = c("all_cause_deaths", "population"), by=c("adm2_id", "adm2_name")] #total number of deaths by muni
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
gamO3.bas <- gam(all_cause_deaths ~ s(temp_mean_2day, bs= "cs", k = 3) + s(o3_mean_2day, bs= "cs", k = 4) + offset(log(population)),
             family=poisson, data=dt)
plot(gamO3.bas, select = 2, shade = T, rug = T) # residuals = T

#move to adding trend variables- seasonal variables significant and in expected direction (summer and spring with lowest mortality)
  #heat effects look to be strongest. Why? Added random effect on municipality
gamO3.bas2 <- gam(all_cause_deaths ~ s(temp_mean_2day, bs= "cs", k = 10) + s(o3_mean_2day, bs= "cs", k = 10) + as.factor(season) +
                  s(trend,bs="cs", k=67, fx=T) + s(adm2_id, bs = "re") + offset(log(population)),
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
#make the sine variable beforehand


#Make ozone and temperature a tensor product-try month trend instead of season**************************BEST MODEL************************
gamO3.ten.mon <- gam(all_cause_deaths ~  te(temp_mean_2day, o3_mean_2day, k=c(4,5), fx = T) + as.factor(month) +
                   s(trend,bs="cs",k=67,fx=T) + s(adm2_id, bs = "re") + offset(log(population)),
                 family=poisson, data=dt)
plot(gamO3.ten.mon, select = 1, pers = T,ticktype="detailed", theta = -35, too.far = .15, se = T, col = "orange2",
       main = "Log relative risk of mortality", xlab = "2-days average temperature [°C]", ylab = "2-day average ozone [µg/m³]") # residuals = T
plot(gamO3.ten.mon, all.terms = T)

gamO3.ten.mon2 <- gam(all_cause_deaths ~  te(temp_mean_day, o3_mean_day, k=c(4,5), fx = T) + as.factor(month) +
                       s(trend,bs="cs",k=35,fx=T) + s(adm2_id, bs = "re") + offset(log(population)),
                     family=poisson, data=dt)
plot(gamO3.ten.mon2, select = 1, pers = T,ticktype="detailed", theta = -35, too.far = .15, se = T, col = "orange2",
     main = "Log relative risk of mortality", xlab = "2-days average temperature [°C]", ylab = "2-day average ozone [µg/m³]") # residuals = T
plot(gamO3.ten.mon2, all.terms = T)

################################################BACK TO BASICS#############################################################################
####BASIC MODEL-Showing individual relationships with 2 day mean
gam.all.bas <- gam(all_cause_deaths ~ s(temp_mean_2day,bs="cs") + s(o3_mean_2day,bs="cs") + s(pm25_mean_2day,bs="cs") + as.factor(season) +
                     s(adm2_id, bs = "re") + s(trend2,bs="cs",k=72,fx=T) + offset(log(population)), #as.factor(season), s(trend2,bs="cs",k=75,fx=T) + as.factor(year_id)
                   family=quasipoisson, data=dt)
plot(gam.all.bas, select = 1, xlab = "2 Day Average Temperature (°C)", ylab = "Log RR Mortality", ylim=c(-.1, .29), rug=T)
plot(gam.all.bas, select = 2, xlab = "2 Day Average Ozone Concentration (ppb)", ylab = "Log RR Mortality", ylim=c(-.1, .29), rug=T)
plot(gam.all.bas, select = 3, xlab = "2 Day Average PM 2.5 Concentration (ppb)", ylab = "Log RR Mortality", ylim=c(-.1, .5), rug=T)

##-Showing individual relationships with 14 day mean
gam.all.bas14 <- gam(all_cause_deaths ~ s(temp_mean_14day,bs="cs") + s(o3_mean_14day,bs="cs") + s(pm25_mean_14day,bs="cs") + as.factor(season) +
                     s(adm2_id, bs = "re") + s(trend2,bs="cs",k=72,fx=T) + offset(log(population)), #as.factor(season), s(trend2,bs="cs",k=75,fx=T) + as.factor(year_id)
                   family=quasipoisson, data=dt)
plot(gam.all.bas14, select = 1, xlab = "14 Day Average Temperature (°C)", ylab = "Log RR Mortality", ylim=c(-.15, .29), rug=T)
plot(gam.all.bas14, select = 2, xlab = "14 Day Average Ozone Concentration (ppb)", ylab = "Log RR Mortality", ylim=c(-.2, .29), rug=T)
plot(gam.all.bas14, select = 3, xlab = "14 Day Average PM 2.5 Concentration (ppb)", ylab = "Log RR Mortality", ylim=c(-.25, .25), rug=T)

####################################TENSOR PRODUCTS######################################################################################
#all cause mortality-2 day average
gamO3.ten.simp <- gam(all_cause_deaths ~  te(o3_mean_2day,temp_mean_2day, k=c(3,3), fx = T) + as.factor(season) + s(trend2,k=72,fx=T) +
                        s(adm2_id, bs = "re") + offset(log(population)), # + pm25_mean_day, k=c(5,5), fx = T
                      family=quasipoisson, data=dt)
plot(gamO3.ten.simp, select = 1, pers = T, ticktype="detailed", theta = -35, too.far = .07, se = T, col = "lightblue",
     main = "Log relative risk of mortality", xlab = "2-day Average temperature [°C]", ylab = "2-day Average ozone [µg/m³]") # residuals = T
plot(gamO3.ten.simp, all.terms = T)

gampm25.ten.simp <- gam(all_cause_deaths ~  te(temp_mean_2day, pm25_mean_2day, k=c(3,3), fx = T) + as.factor(season) + s(trend2,k=72,fx=T) +
                        s(adm2_id, bs = "re") + offset(log(population)),
                      family=quasipoisson, data=dt)
plot(gampm25.ten.simp, select = 1, pers = T, ticktype="detailed", theta = 35, too.far = .07, se = T, col = "orange2",
     main = "Log relative risk of mortality", xlab = "2-day Average temperature [°C]", ylab = "2-day Average PM 2.5 [µg/m³]") # residuals = T
plot(gamO3.ten.simp, all.terms = T)

###all cause mortality-14 day average####################
gamO3.ten.simp <- gam(all_cause_deaths ~  te(temp_mean_14day, o3_mean_14day, k=c(3,3), fx = T) + as.factor(season) + s(trend2,k=72,fx=T) +
                        s(adm2_id, bs = "re") + offset(log(population)), # + pm25_mean_day, k=c(5,5), fx = T
                      family=quasipoisson, data=dt)
plot(gamO3.ten.simp, select = 1, pers = T, ticktype="detailed", theta = -45, too.far = .07, se = T, col = "lightblue",
     main = "Log relative risk of mortality", xlab = "14-day Average temperature [°C]", ylab = "14-day Average ozone [µg/m³]") # residuals = T

gampm25.ten.simp <- gam(all_cause_deaths ~  te(temp_mean_14day, pm25_mean_14day) + as.factor(season) + s(trend2,k=72,fx=T) +
                          s(adm2_id, bs = "re") + offset(log(population)),
                        family=quasipoisson, data=dt)
plot(gampm25.ten.simp, select = 1, pers = T, ticktype="detailed", theta = -25, too.far = .07, se = T, col = "orange2",
     main = "Log relative risk of mortality", xlab = "14-day Average temperature [°C]", ylab = "14-day Average PM 2.5 [µg/m³]") # residuals = T

###cardiovascular cause mortality-2 day average##############
gamO3.ten.simp <- gam(cardio_deaths ~  te(temp_mean_2day, o3_mean_2day, k=c(3,3), fx = T) + as.factor(season) + s(trend2,k=72,fx=T) +
                        s(adm2_id, bs = "re") + offset(log(population)), # + pm25_mean_day, k=c(5,5), fx = T
                      family=quasipoisson, data=dt)
plot(gamO3.ten.simp, select = 1, pers = T, ticktype="detailed", theta = -35, too.far = .07, se = T, col = "lightblue",
     main = "Log relative risk of CVD mortality", xlab = "2-day Average temperature [°C]", ylab = "2-day Average ozone [µg/m³]") # residuals = T

gampm25.ten.simp <- gam(cardio_deaths ~  te(temp_mean_2day, pm25_mean_2day, k=c(3,3), fx = T) + as.factor(season) + s(trend2,k=72,fx=T) +
                          s(adm2_id, bs = "re") + offset(log(population)),
                        family=quasipoisson, data=dt)
plot(gampm25.ten.simp, select = 1, pers = T, ticktype="detailed", theta = -35, too.far = .07, se = T, col = "orange2",
     main = "Log relative risk of CVD mortality", xlab = "2-day Average temperature [°C]", ylab = "2-day Average PM 2.5 [µg/m³]") # residuals = T

####cardiovascular cause mortality-14 day average############
gamO3.ten.simp <- gam(cardio_deaths ~  te(temp_mean_14day, o3_mean_14day, k=c(3,3), fx = T) + as.factor(season) + s(trend2,k=72,fx=T) +
                        s(adm2_id, bs = "re") + offset(log(population)), # + pm25_mean_day, k=c(5,5), fx = T
                      family=quasipoisson, data=dt)
plot(gamO3.ten.simp, select = 1, pers = T, ticktype="detailed", theta = -35, too.far = .07, se = T, col = "lightblue",
     main = "Log relative risk of CVD mortality", xlab = "14-day Average temperature [°C]", ylab = "14-day Average ozone [µg/m³]") # residuals = T

gampm25.ten.simp <- gam(cardio_deaths ~  te(temp_mean_14day, pm25_mean_14day, k=c(3,3), fx = T) + as.factor(season) + s(trend2,k=72,fx=T) +
                          s(adm2_id, bs = "re") + offset(log(population)),
                        family=quasipoisson, data=dt)
plot(gampm25.ten.simp, select = 1, pers = T, ticktype="detailed", theta = -25, too.far = .07, se = T, col = "orange2",
     main = "Log relative risk of CVD mortality", xlab = "14-day Average temperature [°C]", ylab = "14-day Average PM 2.5 [µg/m³]") # residuals = T

#####respiratory cause mortality-2 day average###############
gamO3.ten.simp <- gam(resp_deaths ~  te(temp_mean_2day, o3_mean_2day, k=c(3,3), fx = T) + as.factor(season) + s(trend2,k=72,fx=T) +
                        s(adm2_id, bs = "re") + offset(log(population)), # + pm25_mean_day, k=c(5,5), fx = T
                      family=quasipoisson, data=dt)
plot(gamO3.ten.simp, select = 1, pers = T, ticktype="detailed", theta = -35, too.far = .07, se = T, col = "lightblue",
     main = "Log relative risk of respiratory mortality", xlab = "2-day Average temperature [°C]", ylab = "2-day Average ozone [µg/m³]") # residuals = T

gampm25.ten.simp <- gam(resp_deaths ~  te(temp_mean_2day, pm25_mean_2day, k=c(3,3), fx = T) + as.factor(season) + s(trend2,k=72,fx=T) +
                          s(adm2_id, bs = "re") + offset(log(population)),
                        family=quasipoisson, data=dt)
plot(gampm25.ten.simp, select = 1, pers = T, ticktype="detailed", theta = -35, too.far = .07, se = T, col = "orange2",
     main = "Log relative risk of respiratory mortality", xlab = "2-day Average temperature [°C]", ylab = "2-day Average PM 2.5 [µg/m³]") # residuals = T

######respiratory cause mortality-14 day average################
gamO3.ten.simp <- gam(resp_deaths ~  te(temp_mean_14day, o3_mean_2day, k=c(3,3), fx = T) + as.factor(season) + s(trend2,k=72,fx=T) +
                        s(adm2_id, bs = "re") + offset(log(population)), # + pm25_mean_day, k=c(5,5), fx = T
                      family=quasipoisson, data=dt)
plot(gamO3.ten.simp, select = 1, pers = T, ticktype="detailed", theta = -25, too.far = .07, se = T, col = "lightblue",
     main = "Log relative risk of respiratory mortality", xlab = "14-day Average temperature [°C]", ylab = "2-day Average ozone [µg/m³]") # residuals = T

gampm25.ten.simp <- gam(resp_deaths ~  te(temp_mean_14day, pm25_mean_14day, k=c(3,3), fx = T) + as.factor(season) + s(trend2,k=72,fx=T) +
                          s(adm2_id, bs = "re") + offset(log(population)),
                        family=quasipoisson, data=dt)
plot(gampm25.ten.simp, select = 1, pers = T, ticktype="detailed", theta = -25, too.far = .07, se = T, col = "orange2",
     main = "Log relative risk of respiratory mortality", xlab = "14-day Average temperature [°C]", ylab = "14-day Average PM 2.5 [µg/m³]") # residuals = T


################################################CAUSE-SPECIFIC#############################################################################
####Showing individual relationships with 2 day mean for respiratory causes
gam.all.bas <- gam(resp_deaths ~ s(temp_mean_2day,bs="cs") + s(o3_mean_2day,bs="cs") + s(pm25_mean_2day,bs="cs") + as.factor(season) +
                     s(adm2_id, bs = "re") + s(trend2,bs="cs",k=72,fx=T) + offset(log(population)), family=quasipoisson, data=dt)
plot(gam.all.bas, select = 1, xlab = "2 Day Average Temperature (°C)", ylab = "Log RR Respiratory Mortality", ylim=c(-.05, .1), rug=T)
plot(gam.all.bas, select = 2, xlab = "2 Day Average Ozone Concentration (ppb)", ylab = "Log RR Respiratory Mortality", ylim=c(-.15, .23), rug=T)
plot(gam.all.bas, select = 3, xlab = "2 Day Average PM 2.5 Concentration (ppb)", ylab = "Log RR Respiratory Mortality", ylim=c(-.1, .15), rug=T)

##-Showing individual relationships with 14 day mean for respiratory causes
gam.all.bas14 <- gam(resp_deaths ~ s(temp_mean_14day,bs="cs") + s(o3_mean_14day,bs="cs") + s(pm25_mean_14day,bs="cs") + as.factor(season) +
                       s(adm2_id, bs = "re") + s(trend2,bs="cs",k=72,fx=T) + offset(log(population)), family=quasipoisson, data=dt)
plot(gam.all.bas14, select = 1, xlab = "14 Day Average Temperature (°C)", ylab = "Log RR Respiratory Mortality", ylim=c(-.35, .29), rug=T)
plot(gam.all.bas14, select = 2, xlab = "14 Day Average Ozone Concentration (ppb)", ylab = "Log RR Respiratory Mortality", ylim=c(-.4, .4), rug=T)
plot(gam.all.bas14, select = 3, xlab = "14 Day Average PM 2.5 Concentration (ppb)", ylab = "Log RR Respiratory Mortality", ylim=c(-.35, .7), rug=T)

####Showing individual relationships with 2 day mean for cardiovascular causes
gam.all.bas <- gam(cardio_deaths ~ s(temp_mean_2day,bs="cs") + s(o3_mean_2day,bs="cs") + s(pm25_mean_2day,bs="cs") + as.factor(season) +
                     s(adm2_id, bs = "re") + s(trend2,bs="cs",k=72,fx=T) + offset(log(population)), family=quasipoisson, data=dt)
plot(gam.all.bas, select = 1, xlab = "2 Day Average Temperature (°C)", ylab = "Log RR Cardiovascular Mortality", ylim=c(-.05, .1), rug=T)
plot(gam.all.bas, select = 2, xlab = "2 Day Average Ozone Concentration (ppb)", ylab = "Log RR Cardiovascular Mortality", ylim=c(-.1, .29), rug=T)
plot(gam.all.bas, select = 3, xlab = "2 Day Average PM 2.5 Concentration (ppb)", ylab = "Log RR Cardiovascular Mortality", ylim=c(-.1, .7), rug=T)

##-Showing individual relationships with 14 day mean for cardiovascular causes
gam.all.bas14 <- gam(cardio_deaths ~ s(temp_mean_14day,bs="cs") + s(o3_mean_14day,bs="cs") + s(pm25_mean_14day,bs="cs") + as.factor(season) +
                       s(adm2_id, bs = "re") + s(trend2,bs="cs",k=72,fx=T) + offset(log(population)), family=quasipoisson, data=dt)
plot(gam.all.bas14, select = 1, xlab = "14 Day Average Temperature (°C)", ylab = "Log RR Cardiovascular Mortality", ylim=c(-.15, .29), rug=T)
plot(gam.all.bas14, select = 2, xlab = "14 Day Average Ozone Concentration (ppb)", ylab = "Log RR Cardiovascular Mortality", ylim=c(-.05, .1), rug=T)
plot(gam.all.bas14, select = 3, xlab = "14 Day Average PM 2.5 Concentration (ppb)", ylab = "Log RR Cardiovascular Mortality", ylim=c(-.05, .1), rug=T)


###############################model interactions-temperature as categorical#####################################################################################
#Generate categories for ozone, temperature and PM 2.5
dt[, o3_cat := cut(o3_mean_day, quantile(o3_mean_day, na.rm = T), labels = c("low", "low-mid", "hi-mid", "hi"))]
dt[, temp_cat := cut(temp_range_day, quantile(temp_mean_day, na.rm = T), labels = c("Low", "Mid-low", "Mid-high", "High"))]
#dt[, temp_cat := cut(temp_mean_day, quantile(temp_mean_2day, probs = c(0,.33,.67, 1), na.rm = T), labels = c("low", "mid", "hi"))]
#dt[, temp_cat := cut(temp_mean_day, quantile(temp_mean_day, probs = c(0,.2,.4,.6,.8,1), na.rm = T), labels = c("a", "b", "c", "d", "e"))]
dt[, pm25_cat := cut(pm25_mean_day, quantile(pm25_mean_day, na.rm = T), labels = c("low", "low-mid", "hi-mid", "hi"))]

#model
gamO3.inter.mon <- gam(all_cause_deaths ~ as.factor(temp_cat) + temp_cat:o3_mean_14day + as.factor(season) +
                          s(trend2,bs="cs",k=72,fx=T) + s(adm2_id, bs = "re") + offset(log(population)),
                       family=quasipoisson, data=dt)
#plot(gamO3.inter.mon, all.terms = T)
#summary(gamO3.inter.mon)

#############################
gampm25.inter.mon <- gam(all_cause_deaths ~ as.factor(temp_cat) + as.factor(temp_cat):pm25_mean_14day + as.factor(season) +
                           s(trend2,bs="cs",k=72,fx=T) + s(adm2_id, bs = "re") + offset(log(population)),
                         family=quasipoisson, data=dt)
#plot(gampm25.inter.mon, all.terms = T)
#summary(gampm25.inter.mon)

# gampm25.inter.mon2 <- gam(all_cause_deaths ~ temp_cat + pm25_mean_14day + temp_cat:pm25_mean_14day + as.factor(month) +
#                            s(adm2_id, bs = "re") + offset(log(population)),
#                          family=quasipoisson, data=dt)
# plot(gampm25.inter.mon2, all.terms = T)


##Plot the coefficients in a pretty way###################################
model_table <- summary(gampm25.inter.mon)
n <- grep("pm25", names(gampm25.inter.mon$coef), value = T)
dt.pm25 <- data.table(variable = c("1", "2", "3", "4", "5"),#c("Low", "Mid-low", "Mid-high", "High"),
                      coef = gampm25.inter.mon$coef[names(gampm25.inter.mon$coef) %in% n] %>% as.vector,
                      Pollutant = "PM 2.5",
                      se = model_table$se[names(model_table$se) %in% n] %>% as.vector)
n <- grep("o3", names(gamO3.inter.mon$coef), value = T)
model_table <- summary(gamO3.inter.mon)
dt.o3 <- data.table(variable = c("1", "2", "3", "4", "5"),#c("Low", "Mid-low", "Mid-high", "High"),
                    coef = gamO3.inter.mon$coef[names(gamO3.inter.mon$coef) %in% n] %>% as.vector,
                    Pollutant = "Ozone",
                    se = model_table$se[names(model_table$se) %in% n] %>% as.vector)

dt.plot <- rbind(dt.pm25, dt.o3)
dt.plot[, lower := exp(coef - (1.96*se))]
dt.plot[, upper := exp(coef + (1.96*se))]
dt.plot[, coef := exp(coef)]
dt.plot[, lower := (lower - 1) * 100]
dt.plot[, upper := (upper - 1) * 100]
dt.plot[, coef := (coef - 1) * 100]

pd <- position_dodge(width=0.2)
ggplot(dt.plot, aes(variable, coef, color = Pollutant)) +
  geom_point(position=pd) +
  geom_errorbar(aes(ymin=lower,ymax=upper), position = pd) +
  theme_bw() + scale_x_discrete(limits=c("1", "2", "3", "4", "5")) +
  ylab("% Increase in Mortality per 1 ppb increase in Pollutant") + xlab("Temperature Quantiles") +
  scale_color_manual(values=c("#56B4E9", "orange")) + theme(text = element_text(size=20))
#########################################################################

##model interactions-pollutants as categorical###################################################################################
#create breakpoint variables
bp.heat <- quantile(dt$temp_mean_day, probs = .95, na.rm = T)
bp.cold <- quantile(dt$temp_mean_day, probs = .50, na.rm = T)
dt[, temp01_heat := ifelse(temp_mean_day <= bp.heat, 0, temp_mean_day - bp.heat)]
dt[, temp01_cold := ifelse(temp_mean_day <= bp.cold, 0, temp_mean_day - bp.cold)]

#heat effects
gam.heat.o3.inter <- gam(all_cause_deaths ~ as.factor(o3_cat) + temp_mean_day + temp01_heat:as.factor(o3_cat) + as.factor(season) +
                           s(trend2,bs="cs",k=79,fx=T) + s(adm2_id, bs = "re") + offset(log(population)),
                         family=quasipoisson, data=dt)
plot(gam.heat.o3.inter, all.terms = T)

gam.heat.pm25.inter <- gam(all_cause_deaths ~ as.factor(pm25_cat) + temp_mean_day + temp01_heat:as.factor(pm25_cat) + as.factor(season) +
                           s(trend2,bs="cs",k=79,fx=T) + s(adm2_id, bs = "re") + offset(log(population)),
                         family=quasipoisson, data=dt)
plot(gam.heat.pm25.inter, all.terms = T)
summary(gam.heat.pm25.inter)
##show the cutoffs for PM 2.5 quantiles

#cold effects
gam.cold.o3.inter <- gam(all_cause_deaths ~ as.factor(o3_cat) + temp_mean_day + temp01_cold:as.factor(o3_cat) + as.factor(season) +
                           s(trend2,bs="cs",k=79,fx=T) + s(adm2_id, bs = "re") + offset(log(population)),
                         family=quasipoisson, data=dt)
plot(gam.cold.o3.inter, all.terms = T)

gam.cold.pm25.inter <- gam(all_cause_deaths ~ as.factor(pm25_cat) + temp01_cold + temp_mean_day:as.factor(pm25_cat) + as.factor(season) +
                           s(trend2,bs="cs",k=79,fx=T) + s(adm2_id, bs = "re") + offset(log(population)),
                         family=quasipoisson, data=dt)
summary(gam.cold.pm25.inter)
(as.numeric(exp(summary(gam.cold.pm25.inter)$p.coeff[c(9:12)]))-1) * 100

##Plot the coefficients in a pretty way###################################
model_table <- summary(gam.heat.pm25.inter)
n <- grep("temp01", names(gam.heat.pm25.inter$coef), value = T)
dt.heat <- data.table(variable = c("Low", "Mid-low", "Mid-high", "High"),
                      coef = gam.heat.pm25.inter$coef[names(gam.heat.pm25.inter$coef) %in% n] %>% as.vector,
                      Effect = "Heat",
                      se = model_table$se[names(model_table$se) %in% n] %>% as.vector)
n <- grep("temp_mean", names(gam.cold.pm25.inter$coef), value = T)
model_table <- summary(gam.cold.pm25.inter)
dt.cold <- data.table(variable = c("Low", "Mid-low", "Mid-high", "High"),
                    coef = gam.cold.pm25.inter$coef[names(gam.cold.pm25.inter$coef) %in% n] %>% as.vector,
                    Effect = "Cold",
                    se = model_table$se[names(model_table$se) %in% n] %>% as.vector)

dt.plot <- rbind(dt.heat, dt.cold)
dt.plot[, lower := exp(coef - (1.96*se))]
dt.plot[, upper := exp(coef + (1.96*se))]
dt.plot[, coef := exp(coef)]

pd <- position_dodge(width=0.2)
ggplot(dt.plot, aes(variable, coef, color = Effect)) +
  geom_point(position=pd) +
  geom_errorbar(aes(ymin=lower,ymax=upper), position = pd) +
  theme_bw() + scale_x_discrete(limits=c("Low", "Mid-low", "Mid-high", "High")) +
  ylab("% Increase in Mortality per 1 °C decrease in Temperature") + xlab("PM 2.5 Quantiles")

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

cod <- sql_query(dbname="ghdx",
                 host="ghdx-db-pi01.ihme.washington.edu",
                 query=paste0("SELECT gn.nid, gn.title, CONCAT( 'http://internal-ghdx.healthdata.org/node/', gn.nid ) AS URL, gfdfms.field_microdata_status_tid,gttd.name AS 'Microdata Status' ",
                              "FROM ghdx.node gn ",
                              "LEFT JOIN ghdx.field_data_field_microdata_status gfdfms ON gn.nid = gfdfms.entity_id ",
                              "LEFT JOIN ghdx.taxonomy_term_data gttd ON gfdfms.field_microdata_status_tid = gttd.tid ",
                              "WHERE gn.type = 'record' AND name = 'Have microdata'"))
