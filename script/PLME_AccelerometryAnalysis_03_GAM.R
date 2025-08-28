## Install package with specific version
# remotes::install_version("mgcv", version = "1.9-1") ## require version from 1.9-1

## Load library ----
library(tidyverse)
library(mgcv) 
library(DHARMa)
library(gratia)
library(boot)

## Load data ----
#### Biometric data
meta <- read.csv("data/PLME_accelerometry_meta.csv", header = TRUE, na = c("", "NA"))
meta <- meta[,names(meta) %in% c("turtle.ID", "acce.roundNo", "acce.netWeight",
                                 "sex", "cl", "pl", "wt", "hw", "tl", "leech")]

#### Activity profile data
accel <- read.csv("data/PLME_accelerometry_hmm.csv", header = TRUE)
accel$date <- as.Date(accel$date, format = "%Y-%m-%d")
accel$time <- as.POSIXct(accel$time, format = "%Y-%m-%d %H:%M:%S")
accel$state <- as.factor(accel$state)
str(accel)

####accel#### Temperature data
temp <- read.csv("data/temp_kfbg_daily.csv", header = TRUE)
temp$time <- as.Date(temp$time, format = "%Y-%m-%d")
str(temp)

#### Precipitation data 
prec <- read.csv("data/rainfall_kfbg_daily.csv", header = TRUE)
prec$date <- as.Date(prec$date, format = "%Y-%m-%d")
str(prec)

## Data processing ---- 

### Proportion of time being stationary
accel %>% 
  bind_rows() %>% 
  group_by(turtle.ID, state) %>% 
  summarise(n = n()) %>% 
  complete(state, fill = list(n = 0)) %>%
  mutate(prop = n / sum(n)) %>% 
  filter(state == "1") %>% 
  as.data.frame() %>% 
  summarise(mean = mean(prop), 
            sd = sd(prop))

### Converting to proportion data per hour
state.per.hour <- accel %>% 
  bind_rows() %>% 
  group_by(ID, hour, state) %>% 
  summarise(n = n()) %>% 
  complete(state, fill = list(n = 0)) %>%
  mutate(prop = n / sum(n))  
state.per.hour$turtle.ID <- gsub("(.+).*-(.+).*-(.+)", "\\1", state.per.hour$ID)
state.per.hour$acce.roundNo <- gsub("(.+).*-(.+).*-(.+)", "\\2", state.per.hour$ID)
state.per.hour$track.ID <- gsub("(.+).*-(.+).*-(.+)", "\\3", state.per.hour$ID)

state.per.hour <- merge(meta, state.per.hour, by = c("turtle.ID", "acce.roundNo"))
state.per.hour$turtle.ID <- as.factor(state.per.hour$turtle.ID)
state.per.hour$acce.roundNo <- as.factor(state.per.hour$acce.roundNo)
state.per.hour$track.ID <- as.factor(state.per.hour$track.ID)
state.per.hour$sex <- as.factor(state.per.hour$sex)
state.per.hour$cl <- scale(state.per.hour$cl)
state.per.hour$wt <- scale(state.per.hour$wt)
state.per.hour$acce.netWeight <- scale(state.per.hour$acce.netWeight)
str(state.per.hour)

### Converting to proportion data per day
state.per.day <- accel %>% 
  bind_rows() %>% 
  group_by(ID, date, state) %>% 
  summarise(n = n()) %>% 
  complete(state, fill = list(n = 0)) %>%
  mutate(prop = n / sum(n))
state.per.day$turtle.ID <- gsub("(.+).*-(.+).*-(.+)", "\\1", state.per.day$ID)
state.per.day$acce.roundNo <- gsub("(.+).*-(.+).*-(.+)", "\\2", state.per.day$ID)
state.per.day$track.ID <- gsub("(.+).*-(.+).*-(.+)", "\\3", state.per.day$ID)

state.per.day$julian <- format(state.per.day$date, "%j")
state.per.day$julian <- as.numeric(state.per.day$julian)

state.per.day <- merge(meta, state.per.day, by = c("turtle.ID", "acce.roundNo"))
state.per.day <- merge(state.per.day, temp, by.x = "date", by.y = "time")
state.per.day <- merge(state.per.day, prec, by = "date")
state.per.day$turtle.ID <- as.factor(state.per.day$turtle.ID)
state.per.day$acce.roundNo <- as.factor(state.per.day$acce.roundNo)
state.per.day$track.ID <- as.factor(state.per.day$track.ID)
state.per.day$sex <- as.factor(state.per.day$sex)
state.per.day$cl <- scale(state.per.day$cl)
state.per.day$wt <- scale(state.per.day$wt)
state.per.day$airTemp <- scale(state.per.day$airTemp)
state.per.day$waterTemp <- scale(state.per.day$waterTemp)
state.per.day$rainfall <- scale(state.per.day$rainfall)
state.per.day$acce.netWeight <- scale(state.per.day$acce.netWeight)
str(state.per.day)

## GAM ----

### Diurnal Changes ----

#### Fitting Model
gam.hour <- gam(prop ~ sex + as.vector(wt) +  
                       s(hour, bs = "cc", by = factor(sex), k = -1) + 
                       s(turtle.ID, bs = "re", k = -1) + 
                       s(acce.roundNo, bs = "re", k = -1) +     
                       s(turtle.ID, acce.roundNo, bs = "re", k = -1),  
                data = data.frame(filter(state.per.hour, state == "2")), 
                method = "REML", family = tw(link="logit"), 
                knots = list(hour = c(-0.5, 23.5)))
summary(gam.hour)

#### Model evaluation
gam.check(gam.hour)

gam.hour.simRes <- simulateResiduals(gam.hour)
plot(gam.hour.simRes)

acf(resid(gam.hour))

par(mfrow = c(2,2))
plot(gam.hour, se = TRUE)

#### Plotting -- Partial Effect of Hours
#### Fig 4
smooth_estimates(gam.hour) %>%
  add_confint() %>% 
  filter(.type == "Cyclic CRS") %>% 
  ggplot(aes(y = (.estimate), x = hour)) + 
  geom_ribbon(aes(ymin = (.lower_ci), ymax = (.upper_ci), 
                  col = `factor(sex)`, fill = `factor(sex)`), 
              alpha = 0.05, lty = 2) + 
  geom_line(aes(col = `factor(sex)`)) + 
  scale_color_manual(labels = c("Female", "Male"), values = c("#f8766d", "#00b0f6")) +
  scale_fill_manual(labels = c("Female", "Male"), values = c("#f8766d", "#00b0f6")) +
  ylab("Partial effect on time activity budget") + 
  xlab("Hour") + 
  theme_classic() + 
  theme(legend.title = element_blank(), 
        legend.background = element_rect(fill = "transparent", colour = NA),
        plot.background = element_rect(fill = "transparent", colour = NA), 
        panel.background = element_rect(fill = "transparent", colour = NA)) + 
  theme(text = element_text(size = 12), 
        axis.text = element_text(size = 12), 
        legend.text = element_text(size = 12))
  
#### Daily Changes ---- 

#### Fitting Model 
gam.day <- gam(prop ~ sex + as.vector(wt) +
                      as.vector(airTemp) + as.vector(rainfall) + 
                      s(julian, bs = "cs", k = -1, by = factor(sex)) +
                      s(turtle.ID, bs = "re", k = -1) + 
                      s(acce.roundNo, bs = "re", k = -1) + 
                      s(turtle.ID, acce.roundNo, bs = "re", k = -1),
               data = filter(state.per.day, state == "2"), 
               method = "REML", family = tw(link = "logit"))
summary(gam.day)

#### Model evaluation
gam.check(gam.day)

gam.day.simRes <- simulateResiduals(gam.day)
plot(gam.day.simRes)

acf(resid(gam.day))

par(mfrow = c(2,2))
plot(gam.day, se = TRUE)

#### Prediction ---- 
#### Create new dataset for prediction 
gam.day.newData <- state.per.day %>% filter(state == "2")
gam.day.airTemp <- with(gam.day.newData, 
                        expand.grid(airTemp = seq(min(airTemp), max(airTemp), 0.1),
                                    rainfall = 0, 
                                    sex = factor(levels(sex)), 
                                    wt = 0))
gam.day.rainfall <- with(gam.day.newData, 
                         expand.grid(airTemp = 0,
                                     rainfall = seq(min(rainfall), max(rainfall), 0.1), 
                                     sex = factor(levels(sex)), 
                                     wt = 0))

#### Make prediction
gam.day.airTempPred <- predict(gam.day, gam.day.airTemp, 
                               exclude = c('s(julian):factor(sex)F',
                                           's(julian):factor(sex)M',
                                           's(turtle.ID)', 
                                           's(acce.roundNo)', 
                                           's(turtle.ID,acce.roundNo)'), 
                               newdata.guaranteed = TRUE, 
                               se = TRUE)
gam.day.airTemp <- cbind(gam.day.airTemp, gam.day.airTempPred)

gam.day.airTemp$upper_ci <- gam.day.airTemp$fit + gam.day.airTemp$se.fit
gam.day.airTemp$lower_ci <- gam.day.airTemp$fit - gam.day.airTemp$se.fit

gam.day.airTemp$fit <- inv.logit(gam.day.airTemp$fit)
gam.day.airTemp$upper_ci <- inv.logit(gam.day.airTemp$upper_ci)
gam.day.airTemp$lower_ci <- inv.logit(gam.day.airTemp$lower_ci)

gam.day.rainfallPred <- predict(gam.day, gam.day.rainfall, 
                                exclude = c('s(julian):factor(sex)F',
                                            's(julian):factor(sex)M',
                                            's(turtle.ID)', 
                                            's(acce.roundNo)', 
                                            's(turtle.ID,acce.roundNo)'), 
                                newdata.guaranteed = TRUE, 
                                se = TRUE)
gam.day.rainfall <- cbind(gam.day.rainfall, gam.day.rainfallPred)

gam.day.rainfall$upper_ci <- gam.day.rainfall$fit + gam.day.rainfall$se.fit
gam.day.rainfall$lower_ci <- gam.day.rainfall$fit - gam.day.rainfall$se.fit

gam.day.rainfall$fit <- inv.logit(gam.day.rainfall$fit)
gam.day.rainfall$upper_ci <- inv.logit(gam.day.rainfall$upper_ci)
gam.day.rainfall$lower_ci <- inv.logit(gam.day.rainfall$lower_ci)

#### Plotting ----
#### Fig S3
gam.day.airTemp %>% 
  filter(sex == "F") %>% 
  ggplot(aes(x = airTemp * attr(gam.day.newData$airTemp, 'scaled:scale') + attr(gam.day.newData$airTemp, 'scaled:center'), 
             y = fit)) + 
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), 
              alpha = 0.1, lty = 2, col = "#1cade4", fill = "#1cade4") +
  geom_line(col = "#1cade4") +
  ylab("Proportion of time spent being active") + 
  xlab("Ambient temperature (°C)") + 
  theme_classic() + 
  theme(legend.title = element_blank(), 
        legend.background = element_rect(fill = "transparent", colour = NA), 
        plot.background = element_rect(fill = "transparent", colour = NA), 
        panel.background = element_rect(fill = "transparent", colour = NA)) + 
  theme(text = element_text(size = 12), 
        axis.text = element_text(size = 12))

#### Fig S4
gam.day.rainfall %>% 
  filter(sex == "F") %>% 
  ggplot(aes(x = rainfall * attr(gam.day.newData$rainfall, 'scaled:scale') + attr(gam.day.newData$rainfall, 'scaled:center'), 
             y = fit)) + 
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci),
              alpha = 0.1, lty = 2, col = "#1cade4", fill = "#1cade4") +
  geom_line(col = "#1cade4") +
  ylab("Proportion of time spent being active") + 
  xlab("Daily precipitation (mm)") + 
  theme_classic() + 
  theme(legend.title = element_blank(), 
        legend.background = element_rect(fill = "transparent", colour = NA), 
        plot.background = element_rect(fill = "transparent", colour = NA), 
        panel.background = element_rect(fill = "transparent", colour = NA)) + 
  theme(text = element_text(size = 12), 
        axis.text = element_text(size = 12))
