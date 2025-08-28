# Load packages ----
library(tidyverse)
library(xts)
library(highfrequency)
library(RcppRoll)
library(TSstudio)
library(patchwork)
library(PerformanceAnalytics)
library(usdm)

# Load data ----
## Load metadata ----
meta <- read.csv("data/PLME_accelerometry_meta.csv", header = TRUE, na = c("", "NA"))
meta$acce.startDate <- sub("([A-Za-z]{3})[a-z]", "\\1", meta$acce.startDate)
meta$acce.startTime <- sub("([A-Za-z]{3})[a-z]", "\\1", meta$acce.startTime)
str(meta)

## Load accelerometer data ----
accel_file <- list.files(path = "data/accelerometer", 
                         pattern = "PLME.*.data.*csv",
                         full.names = TRUE, 
                         recursive = TRUE)

accel <- lapply(accel_file, function(x){
  dat <- read.csv(x, header = FALSE, na = c("")) %>% drop_na()
  names(dat) <- c("time", "x", "y", "z", "temp")
  dat$filename <- as.character(x)
  return(dat)
})
accel <- do.call(rbind, accel)
accel$y <- as.integer(accel$y)
str(accel)

accel <- merge(meta, accel, by = c("filename"))
accel$turtle.ID <- as.factor(accel$turtle.ID)
accel$acce.roundNo <- as.factor(accel$acce.roundNo)
accel$filename <- as.factor(accel$filename)

## Convert counting time in milliseconds to datetime format
my_options <- options(digits.secs = 3)

accel$acce.startTime <- as.POSIXct(accel$acce.startTime, "%d/%m/%Y %H:%M", tz = "GMT")
accel$time <- as.character(accel$time)
accel$time <- str_split(accel$time, " : ", simplify = TRUE)
accel$time <- accel$acce.startTime + 
  days(accel$time[,1]) + 
  hours(accel$time[,2]) + 
  minutes(accel$time[,3]) + 
  seconds(accel$time[,4])

## Unit conversion (g)
accel$x <- accel$x*(4/256)
accel$y <- accel$y*(4/256)
accel$z <- accel$z*(4/256)

## Load temperature data ----
temp_file <- list.files(path = "data/thermochrone", 
                        pattern = ".*csv", 
                        full.names = TRUE, 
                        recursive = TRUE)

temp <- lapply(temp_file, function(x){
  dat <- read.csv(x, skip = 14, header = TRUE, row.names = NULL) %>% drop_na()
  names(dat) <- c("date", "time", "unit", "temp")
  dat$filename <- as.character(x)
  dat$medium <- gsub(".*stream_(.+)Temp.*", "\\1", x)
  dat$loc <- gsub(".*KFBG_(.+)stream.*", "\\1", x)
  dat$round <- gsub(".*round(.+)/.*", "\\1", x)
  
  dat$time <- paste(dat$date, dat$time, " ")
  if(nchar(gsub("(.+).*/(.+).*/(.+)", "\\3", dat$date[1])) == 4){
    if(length(unique(gsub("(.+).*/(.+).*/(.+)", "\\2", dat$date)))<=12){
      dat$time <- as.POSIXct(dat$time, format = "%d/%m/%Y %I:%M:%S %p", tz = "GMT")
    }
    else{
      dat$time <- as.POSIXct(dat$time, format = "%m/%d/%Y %I:%M:%S %p", tz = "GMT")
    }
  }
  if(nchar(gsub("(.+).*/(.+).*/(.+)", "\\3", dat$date[1])) == 2){
    if(length(unique(gsub("(.+).*/(.+).*/(.+)", "\\2", dat$date)))<=12){
      dat$time <- as.POSIXct(dat$time, format = "%d/%m/%y %I:%M:%S %p", tz = "GMT")
    }
    else{
      dat$time <- as.POSIXct(dat$time, format = "%m/%d/%y %I:%M:%S %p", tz = "GMT")
    }
  }

  return(dat)
})
temp <- do.call(rbind, temp)

temp$medium <- as.factor(temp$medium)
temp$loc <- as.factor(temp$loc)
temp$round <- as.factor(temp$round)
str(temp)

## Load precipitation data ---- 
rainfall <- read.csv("data/rainfall_kfbg_hourly.csv", header = TRUE, na = "n/a")

rainfall$Date <- as.Date(as.character(rainfall$Date), format = "%Y%m%d")
rainfall$hour <- rainfall$hour - 1
rainfall$Rainfall..mm. <- gsub("[*]", "", rainfall$Rainfall..mm.)
rainfall$Rainfall..mm. <- as.numeric(rainfall$Rainfall..mm.)

rainfall <- rainfall %>% 
  rename(date = Date) %>% 
  rename(rainfall = Rainfall..mm.)

rainfall$roll_sum_24h <- roll_sum(rainfall$rainfall, 24, align = "right", fill = NA)

rainfall$heavyRain <- ifelse(rainfall$roll_sum_24h > 30, "yes", "no")
rainfall$heavyRain <- as.factor(rainfall$heavyRain)

str(rainfall)

## Summarize daily total precipitation
rainfall.day <- rainfall %>% 
  group_by(date) %>% 
  summarise(rainfall = sum(rainfall))
# write.csv(rainfall.day, "data/rainfall_kfbg_daily.csv", row.names = FALSE)

# Time-series analysis ----
## Convert df to xts object ----
#### Acceleration data
accel$turtle.round.ID <- paste0(accel$turtle.ID, "-", accel$acce.roundNo)
accel$turtle.round.ID <- as.factor(accel$turtle.round.ID)
accel_iid <- split(accel, accel$turtle.round.ID)
accelTS_iid <- lapply(accel_iid, function(x){
  as.xts(x[, c("x", "y", "z", "temp")], order.by = x[, c("time")])
})

#### Temperature data
temp$loc.round <- paste(temp$round, temp$loc, temp$medium, sep = "-")
temp$loc.round <- as.factor(temp$loc.round)
temp_rd <- split(temp, temp$loc.round)
tempTS_rd <- lapply(temp_rd, function(x){
  as.xts(x[, c("temp")], order.by = x[, c("time")])
})

## Aggregate time series data ----
aggregateAccel = function(tsData, dur, per, fun = "previoustick"){
  start = min(index(tsData))
  end = max(index(tsData))
  
  if(end - start >= dur){
    
    aggregateDf <- aggregateTS(ts = tsData, alignBy = per, alignPeriod = as.numeric(dur, per), FUN = fun)
    ts_plot(aggregateDf)
    
    return(aggregateDf)
    
  }
  
}
#### Acceleration data
accelTS15min <- lapply(accelTS_iid, aggregateAccel, dur = minutes(15), per = "minutes", fun = "mean")

#### Temperature data
tempTS15min <- lapply(tempTS_rd, aggregateAccel, dur = minutes(15), per = "minutes", fun = "mean")

tempTS1hrs <- lapply(tempTS_rd, aggregateAccel, dur = hours(1), per = "hours", fun = "mean")
tempTS1day <- lapply(tempTS_rd, aggregateAccel, dur = days(1), per = "days", fun = "mean")

## Plot aggregated time series data ----
#### Style 1: Plot of x, y, z in the same panel 
#### Style 2: Plot of x, y, z and temp in individual panels

plotAccel = function(tsName, tsList, style, saved){
  
  if(length(tsList[[tsName]] != 0)){

    start = min(index(tsList[[tsName]]))
    end = max(index(tsList[[tsName]]))
    
    min_y = min(pretty(c(tsList[[tsName]]$x, tsList[[tsName]]$y, tsList[[tsName]]$z)))
    max_y = max(pretty(c(tsList[[tsName]]$x, tsList[[tsName]]$y, tsList[[tsName]]$z)))
    by_y = pretty(c(tsList[[tsName]]$x, tsList[[tsName]]$y, tsList[[tsName]]$z))[2] - pretty(c(tsList[[tsName]]$x, tsList[[tsName]]$y, tsList[[tsName]]$z))[1]
    label = paste0("PLME ", gsub("-.*", "", tsName), " - Rd ", gsub(".*-", "", tsName))
    
    if(style == 1){
      fullAccel <- ggplot(tsList[[tsName]], aes(x = Index, y = x)) + 
        geom_line(col = "#FF0000") +
        geom_line(aes(y = y), col = "#CC6600") +
        geom_line(aes(y = z), col = "#669900") +
        scale_x_datetime(breaks = seq(from = as.POSIXct(round_date(start, "day")), 
                                      to = as.POSIXct(round_date(end, "day")), 
                                      by = "24 hours"), 
                         date_labels = "%b %d %H:%M", 
                         minor_breaks = seq(from = as.POSIXct(round_date(start, "day")), 
                                            to = as.POSIXct(round_date(end, "day")), 
                                            by = "12 hours")) + 
        
        theme_classic() + 
        ylab("Acceleration (g)") + 
        scale_y_continuous(limits = c(min_y, max_y), 
                           breaks = seq(from = min_y, to = max_y, by = by_y)) +
        theme(axis.title.x = element_blank()) +
        theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0)) + 
        geom_text(x = as.POSIXct(round_date(start, "day")), 
                  y = max_y + .1,
                  label = label)
      fullAccel
      
    }
    
    if(style == 2){
      xAccel <- ggplot(tsList[[tsName]], aes(x = Index, y = x)) +
        geom_line(col = "#FF0000") +
        scale_x_datetime(breaks = seq(from = as.POSIXct(round_date(start, "day")), 
                                      to = as.POSIXct(round_date(end, "day")), 
                                      by = "24 hours"), 
                         date_labels = "%b %d %H:%M", 
                         minor_breaks = seq(from = as.POSIXct(round_date(start, "day")), 
                                            to = as.POSIXct(round_date(end, "day")), 
                                            by = "12 hours")) + 
        theme_classic() + 
        ylab("Surge (g)") + 
        scale_y_continuous(limits = c(min_y, max_y), 
                           breaks = seq(from = min_y, to = max_y, by = by_y)) +
        theme(axis.title.x = element_blank()) +
        theme(axis.text.x = element_blank()) + 
        theme(axis.line.x = element_blank()) + 
        theme(axis.ticks.x = element_blank())
      
      yAccel <- ggplot(tsList[[tsName]], aes(x = Index, y = y)) +
        geom_line(col = "#CC6600") +
        scale_x_datetime(breaks = seq(from = as.POSIXct(round_date(start, "day")), 
                                      to = as.POSIXct(round_date(end, "day")), 
                                      by = "24 hours"), 
                         date_labels = "%b %d %H:%M", 
                         minor_breaks = seq(from = as.POSIXct(round_date(start, "day")), 
                                            to = as.POSIXct(round_date(end, "day")), 
                                            by = "12 hours")) + 
        theme_classic() + 
        ylab("Sway (g)") + 
        scale_y_continuous(limits = c(min_y, max_y), 
                           breaks = seq(from = min_y, to = max_y, by = by_y)) +
        theme(axis.title.x = element_blank()) +
        theme(axis.text.x = element_blank()) + 
        theme(axis.line.x = element_blank()) + 
        theme(axis.ticks.x = element_blank())
      
      zAccel <- ggplot(tsList[[tsName]], aes(x = Index, y = z)) +
        geom_line(col = "#669900") +
        scale_x_datetime(breaks = seq(from = as.POSIXct(round_date(start, "day")), 
                                      to = as.POSIXct(round_date(end, "day")), 
                                      by = "24 hours"), 
                         date_labels = "%b %d %H:%M", 
                         minor_breaks = seq(from = as.POSIXct(round_date(start, "day")), 
                                            to = as.POSIXct(round_date(end, "day")), 
                                            by = "12 hours")) + 
        theme_classic() + 
        ylab("Heave (g)") + 
        scale_y_continuous(limits = c(min_y, max_y), 
                           breaks = seq(from = min_y, to = max_y, by = by_y)) +
        theme(axis.title.x = element_blank()) +
        theme(axis.text.x = element_blank()) + 
        theme(axis.line.x = element_blank()) + 
        theme(axis.ticks.x = element_blank())
      
      temp <- ggplot(tsList[[tsName]], aes(x = Index, y = temp)) +
        geom_line() +
        scale_x_datetime(breaks = seq(from = as.POSIXct(round_date(start, "day")), 
                                      to = as.POSIXct(round_date(end, "day")), 
                                      by = "24 hours"), 
                         date_labels = "%b %d %H:%M", 
                         minor_breaks = seq(from = as.POSIXct(round_date(start, "day")), 
                                            to = as.POSIXct(round_date(end, "day")), 
                                            by = "12 hours")) + 
        theme_classic() + 
        ylab("Temperature") + 
        theme(axis.title.x = element_blank()) +
        theme(axis.text.x = element_text(angle = 25, vjust = 1.0, hjust = 1.0))
      
        fullAccel <- xAccel + yAccel + zAccel + temp +
          plot_layout(ncol = 1, heights = c(3,3,3,1)) +
          plot_annotation(title = label)
        fullAccel

    }
    
    if(saved == TRUE){
      
      ggsave(filename = paste0("figures/PLME_", gsub("-.*", "", tsName), "_Rd ", gsub(".*-", "", tsName), "_Accel15min.png"), 
             fullAccel, 
             height = 23.82, width = 24.64, units = "cm")
      
    }
    
    return(fullAccel)
    
  }
  
}

lapply(names(accelTS15min), plotAccel, tsList = accelTS15min, style = 2, saved = FALSE)

## Compute equal interval time-series data ----
computeRawAccel = function(tsData, dur, per){
  
  if(length(tsData) != 0){
    
    start = min(index(tsData))
    end = max(index(tsData))
    
    if(end - start >= dur){
      aggregateDf <- aggregateTS(ts = tsData, alignBy = per, alignPeriod = as.numeric(dur, per), 
                                 FUN = "previoustick", dropna = TRUE) 
      empty <- xts(NULL, seq(start(aggregateDf), end(aggregateDf), by = "1 sec"))
      
      aggregateDf <- merge(aggregateDf, empty, all = TRUE)
      return(aggregateDf)
      
    }
    
  }
  
}

accelTS1sec <- lapply(accelTS_iid, computeRawAccel, dur = seconds(1), per = "seconds")

## Static Body Acceleration ----
staticAccelTS1sec <- lapply(accelTS1sec, function(x){
  rollapply(x, 2, mean, na.rm = TRUE, 
            align = "right", fill = NA) 
})

## Dynamic Body Acceleration ----
dba = function(rawAccelList, staticAccelList){
  
  if(nrow(rawAccelList) == nrow(staticAccelList)){
    
    dynamicAccel <- rawAccelList[, c("x", "y", "z")] - staticAccelList[, c("x", "y", "z")]
    return(dynamicAccel)
    
  }
  else{
    
    print("error: lengths of `rawAccel` and `staticAccel` are not equal.")
    
  }
  
}
dynamicAccelTS1sec <- mapply(dba, rawAccelList = accelTS1sec, staticAccelList = staticAccelTS1sec)

## Overall Dynamic Body Acceleration ----
odba = function(dynamicAccel){
  
  if(length(dynamicAccel) != 0){
    
    dynamicAccel$ODBA <- abs(dynamicAccel$x) + abs(dynamicAccel$y) + abs(dynamicAccel$z)
    return(dynamicAccel)
  }
  
}
dynamicAccelTS1sec <- lapply(dynamicAccelTS1sec, odba)

## Plots of ODBAs ----
staticAccel <- lapply(names(staticAccelTS1sec), 
                      FUN = function(tsName, tsList){
                        
                        tsList[[tsName]]$turtle.ID <- gsub("-.*", "", tsName)
                        tsList[[tsName]]$acce.roundNo <- gsub(".*-", "", tsName)
                        
                        staticAccel <- data.frame(time = index(tsList[[tsName]]), 
                                                  coredata(tsList[[tsName]]))
                        
                        return(staticAccel)
                        
                      }, 
                      tsList = staticAccelTS1sec) %>% bind_rows() %>% drop_na()
dynamicAccel <- lapply(names(dynamicAccelTS1sec), 
                       FUN = function(tsName, tsList){
                         
                         tsList[[tsName]]$turtle.ID <- gsub("-.*", "", tsName)
                         tsList[[tsName]]$acce.roundNo <- gsub(".*-", "", tsName)
                         
                         dynamicAccel <- data.frame(time = index(tsList[[tsName]]), 
                                                    coredata(tsList[[tsName]]))
                         
                         return(dynamicAccel)
                         
                       }, 
                       tsList = dynamicAccelTS1sec) %>% bind_rows() %>% drop_na()
dynamicAccel <- merge(dynamicAccel, staticAccel[, c("time", "turtle.ID", "acce.roundNo")], 
                      by = c("time", "turtle.ID", "acce.roundNo"))

dynamicAccel$turtle.ID <- as.factor(dynamicAccel$turtle.ID)
dynamicAccel$acce.roundNo <- as.factor(dynamicAccel$acce.roundNo)

dynamicAccel <- merge(meta, dynamicAccel, by = c("turtle.ID", "acce.roundNo"))
str(dynamicAccel)

ggplot(dynamicAccel, aes(x = time, y = ODBA)) + 
  theme_bw() +
  theme(panel.grid.major.x = element_line(colour = "black", linewidth = 0.5)) +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.border = element_blank()) +
  theme(panel.spacing = unit(0.2, "lines")) + 
  geom_rect(data = meta[, c("turtle.ID", "sex")] %>% filter(!duplicated(turtle.ID)), aes(fill = sex), 
            xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, alpha = .8, inherit.aes = FALSE) + 
  scale_fill_manual(values = c("#FFE599", "#D8EEFC")) + 
  geom_line(linewidth = 0.8) + 
  facet_grid(turtle.ID ~ acce.roundNo, scales = "free_x", space = "free") +
  scale_x_datetime(breaks = seq(from = as.POSIXct(ceiling_date(min(dynamicAccel$time), "day")), 
                                to = as.POSIXct(floor_date(max(dynamicAccel$time), "day")), 
                                by = "24 hours"), 
                   date_labels = "%b %d %H:%M", 
                   minor_breaks = seq(from = as.POSIXct(ceiling_date(min(dynamicAccel$time), "day")), 
                                      to = as.POSIXct(floor_date(max(dynamicAccel$time), "day")), 
                                      by = "24 hours")) + 
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 40, vjust = 1.0, hjust = 1.0),
        axis.ticks.length.x = unit(.25, "cm"))  + 
  scale_y_continuous(limits = c(min(pretty(dynamicAccel$ODBA)), max(pretty(dynamicAccel$ODBA))),
                     breaks = seq(from = min(pretty(dynamicAccel$ODBA)), to = max(pretty(dynamicAccel$ODBA)), 
                                  by = pretty(dynamicAccel$ODBA)[2] - pretty(dynamicAccel$ODBA)[1])) + 
  theme(axis.text.y = element_blank(), 
        axis.ticks.y = element_blank())  + 
  theme(legend.position = "bottom")

## Data cleaning ----
### Environmental temperature (per 15 mins) ----
temp_mean <- lapply(names(tempTS15min), 
                    FUN = function(tsName, tsList){
                      
                      tempData <- data.frame(time = index(tsList[[tsName]]), 
                                             temp = coredata(tsList[[tsName]]), 
                                             round = gsub("(.+).*-(.+).*-(.+)", "\\1", tsName), 
                                             loc = gsub("(.+).*-(.+).*-(.+)", "\\2", tsName),
                                             medium = gsub("(.+).*-(.+).*-(.+)", "\\3", tsName))
                      
                      return(tempData)
                      
                    }, 
                    tsList = tempTS15min) %>% bind_rows()
airTemp.ups <- temp_mean %>% filter(medium == "Air") %>% filter(loc == "Up")
airTemp.dws <- temp_mean %>% filter(medium == "Air") %>% filter(loc == "Dw")
airTemp <- merge(airTemp.ups[, c("time", "round", "x")], airTemp.dws[, c("time", "round", "x")], 
                 by = c("time", "round"), all = TRUE)
airTemp$airTemp <- apply(airTemp[, c("x.x", "x.y")], 1, mean, na.rm = TRUE)
waterTemp.ups <- temp_mean %>% filter(medium == "Water") %>% filter(loc == "Up")
waterTemp.dws <- temp_mean %>% filter(medium == "Water") %>% filter(loc == "Dw")
waterTemp <- merge(waterTemp.ups[, c("time", "round", "x")], waterTemp.dws[, c("time", "round", "x")], 
                   by = c("time", "round"), all = TRUE)
waterTemp$waterTemp <- apply(waterTemp[, c("x.x", "x.y")], 1, mean, na.rm = TRUE)
temp_mean <- merge(airTemp[, c("time", "round", "airTemp")], waterTemp[, c("time", "round", "waterTemp")], 
                   by = c("time", "round"), all = TRUE)
remove(airTemp.ups, airTemp.dws, airTemp, 
       waterTemp.ups, waterTemp.dws, waterTemp)

### Environmental temperature (per 1 day) ----
dailyTemp_mean <- lapply(names(tempTS1day), 
                         FUN = function(tsName, tsList){
                           
                           tempData <- data.frame(time = index(tsList[[tsName]]), 
                                                  temp = coredata(tsList[[tsName]]), 
                                                  round = gsub("(.+).*-(.+).*-(.+)", "\\1", tsName), 
                                                  loc = gsub("(.+).*-(.+).*-(.+)", "\\2", tsName),
                                                  medium = gsub("(.+).*-(.+).*-(.+)", "\\3", tsName))
                           
                           return(tempData)
                           
                         }, 
                         tsList = tempTS1day) %>% bind_rows()
dailyAirTemp.ups <- dailyTemp_mean %>% filter(medium == "Air") %>% filter(loc == "Up")
dailyAirTemp.dws <- dailyTemp_mean %>% filter(medium == "Air") %>% filter(loc == "Dw")
dailyAirTemp <- merge(dailyAirTemp.ups[, c("time", "round", "x")], dailyAirTemp.dws[, c("time", "round", "x")], 
                      by = c("time", "round"), all = TRUE)
dailyAirTemp$airTemp <- apply(dailyAirTemp[, c("x.x", "x.y")], 1, mean, na.rm = TRUE)
dailyWaterTemp.ups <- dailyTemp_mean %>% filter(medium == "Water") %>% filter(loc == "Up")
dailyWaterTemp.dws <- dailyTemp_mean %>% filter(medium == "Water") %>% filter(loc == "Dw")
dailyWaterTemp <- merge(dailyWaterTemp.ups[, c("time", "round", "x")], dailyWaterTemp.dws[, c("time", "round", "x")], 
                        by = c("time", "round"), all = TRUE)
dailyWaterTemp$waterTemp <- apply(dailyWaterTemp[, c("x.x", "x.y")], 1, mean, na.rm = TRUE)
dailyTemp_mean <- merge(dailyAirTemp[, c("time", "round", "airTemp")], dailyWaterTemp[, c("time", "round", "waterTemp")], 
                        by = c("time", "round"), all = TRUE)
remove(dailyAirTemp.ups, dailyAirTemp.dws, dailyAirTemp, 
       dailyWaterTemp.ups, dailyWaterTemp.dws, dailyWaterTemp)
# write.csv(dailyTemp_mean, "data/temp_kfbg_daily.csv", row.names = FALSE)

### Dynamic Body Acceleration (per 15 mins) ----
dynamicAccel$turtle.round.ID <- paste0(dynamicAccel$turtle.ID, "-", dynamicAccel$acce.roundNo)
dynamicAccel$turtle.round.ID <- as.factor(dynamicAccel$turtle.round.ID)

dynamicAccel_iid <- split(dynamicAccel, dynamicAccel$turtle.round.ID)
dynamicAccelTS15min <- lapply(dynamicAccel_iid, function(x){
  
  #### Convert list of df back into list of xts
  dynamicAccelTS <- as.xts(x[, c("x", "y", "z", "ODBA")], order.by = x[, c("time")])
  
  #### Compute 15-minutes average 
  dynamicAccelTS <- aggregateAccel(dynamicAccelTS, dur = minutes(15), per = "minutes", fun = "mean")
  return(dynamicAccelTS)
  
})

dynamicAccel_mean <- lapply(names(dynamicAccelTS15min), 
                            FUN = function(tsName, tsList){
                              
                              if(length(tsList[[tsName]] != 0)){
                                
                                tsList[[tsName]]$turtle.ID <- gsub("-.*", "", tsName)
                                tsList[[tsName]]$acce.roundNo <- gsub(".*-", "", tsName)
                                
                                dynamicAccel <- data.frame(time = index(tsList[[tsName]]), 
                                                           coredata(tsList[[tsName]]))
                                
                                return(dynamicAccel)
                                
                              }
                              
                            }, 
                            tsList = dynamicAccelTS15min) %>% bind_rows() %>% drop_na(ODBA)

### Combine datasets ----
dynamicAccel_mean$turtle.ID <- as.factor(dynamicAccel_mean$turtle.ID)
dynamicAccel_mean$acce.roundNo <- as.factor(dynamicAccel_mean$acce.roundNo)

dynamicAccel_mean$date <- as.Date(dynamicAccel_mean$time, format = "%Y-%m-%d")
dynamicAccel_mean$year <- format(dynamicAccel_mean$time, "%Y")
dynamicAccel_mean$year <- as.factor(dynamicAccel_mean$year)
dynamicAccel_mean$month <- format(dynamicAccel_mean$time, "%m")
dynamicAccel_mean$month <- as.numeric(dynamicAccel_mean$month)
dynamicAccel_mean$hour <- format(dynamicAccel_mean$time, "%H")
dynamicAccel_mean$hour <- as.numeric(dynamicAccel_mean$hour)
dynamicAccel_mean$julian <- format(dynamicAccel_mean$time, format = "%j")
dynamicAccel_mean$julian <- as.numeric(dynamicAccel_mean$julian)

#### Merge acceleration data with environmental data
dynamicAccel_mean <- merge(dynamicAccel_mean, temp_mean[, !names(temp_mean) %in% c("round")], 
                           by = c("time"), all.x = TRUE)
dynamicAccel_mean <- merge(dynamicAccel_mean, rainfall, 
                           by = c("date", "hour"), all.x = TRUE)

#### Merge acceleration data with biological variables
dynamicAccel_mean <- merge(meta[, !names(meta) %in% c("date")], dynamicAccel_mean, 
                           by = c("turtle.ID", "acce.roundNo"))

str(dynamicAccel_mean)

### Normality check ----
#### Present the histogram
#### Present the Q-Q plot
#### Present the calculated Dmax for K-S test
#### Check the critical (Dmax) alpha,n using D-table
checkNormality <- function(data) {
  print(deparse(substitute(data)))
  plot(density(data))
  hist(data)
  print(paste('Difference between mean and trimmed mean = ',
              mean(data) - mean(data, trim=0.05)))
  print(paste('Standard error = ',sqrt(var(data)/length(data))))
  library(moments)
  print(paste('Skewness = ',skewness(data)))
  print(paste('Kurtosis = ',kurtosis(data)))
  print(ks.test(data,"pnorm", mean=mean(data), sd=sd(data)))
  qqnorm(data); qqline(data)  
}
checkNormality(dynamicAccel_mean$ODBA) ## not normal

### Multicollinearity check ----
#### Continuous variables for checking pairwise correlation
corrCheck <- dynamicAccel_mean %>% 
  dplyr::select(cl, pl, wt, hour, month, airTemp, waterTemp, rainfall) %>% 
  data.frame() %>% 
  drop_na()
corrCheck <- scale(corrCheck)
corrCheck[which(is.na(corrCheck)),] 

#### Check Normal Distribution
apply(corrCheck, 2, checkNormality)

#### Calculate Spearman's Correlation Coefficient
chart.Correlation(corrCheck,
                  method="spearman",
                  histogram=TRUE,
                  pch=16)

#### Calculate VIF
vifstep(corrCheck) ## - pl, - wt, - airTemp

### Saving the prepared data ----

#### Variables to be included 
####
####
#### Biological vars.
#### - sex
#### - cl
####
####
#### Abiotic Vars.
#### - hour
#### - month
#### - season
#### - airTemp / waterTemp
#### - rainfall
####
####
#### Confounding factors
#### - acce.roundNo
#### - acce.netWeight
####
####

md.Data <- dynamicAccel_mean %>% 
  filter(acce.roundNo %in% c("1", "3", "4", "5"))%>%
  dplyr::select(acce.roundNo, 
                turtle.ID, sex, cl, 
                time, date, hour, month, year, julian, 
                x, y, z, ODBA, 
                airTemp, waterTemp, 
                rainfall, roll_sum_24h, heavyRain, 
                acce.netWeight) %>% 
  arrange(time) %>% 
  transform(time = as.character(format(time))) %>% 
  drop_na() 
str(md.Data)

# write.csv(md.Data, "data/PLME_accelerometry_aggregated.csv", row.names = FALSE)


