## Load library ---- 
library(tidyverse)
library(ggpubr)
library(gridExtra)
library(momentuHMM)
library(Rmixmod)
source("script/utility_functions.R")

## Load data ----
md.Data <- read.csv("data/PLME_accelerometry_aggregated.csv", header = TRUE)

md.Data$acce.roundNo <- as.factor(md.Data$acce.roundNo)
md.Data$turtle.ID <- as.factor(md.Data$turtle.ID)
md.Data$ID <- paste(md.Data$turtle.ID, md.Data$acce.roundNo, sep = "-")

md.Data$sex <- as.factor(md.Data$sex)

md.Data$time <- as.POSIXct(md.Data$time, format = "%Y-%m-%d %H:%M:%S")
md.Data$date <- as.Date(md.Data$date, foramt = "%Y-%m-%d")
md.Data$heavyRain <- as.factor(md.Data$heavyRain)
str(md.Data)

## Data Exploration ----
circData <- md.Data %>% 
  group_by(hour, sex) %>% 
  summarise(mean = mean(ODBA))

circData %>% 
  ggplot(aes(x = hour, y = mean, fill = sex)) + 
  facet_grid(~sex, 
             labeller = labeller(sex = c("F" = "Female", 
                                         "M" = "Male"))) +
  geom_bar(stat = "identity") +
  coord_polar(start = 0) + 
  theme_bw() + 
  ylab("Mean ODBA") + 
  scale_fill_manual(values = c("M" = "#084C61", "F" = "#db3a34")) + 
  theme(legend.position = "none") + 
  theme(panel.background = element_rect(fill = "transparent", colour = NA), 
        plot.background = element_rect(fill = "transparent", colour = NA))

## Hidden Markov Models ----
with(md.Data, table(acce.roundNo, turtle.ID))

### Getting the data ready ----
md.Data$airTemp <- scale(md.Data$airTemp)
md.Data$waterTemp <- scale(md.Data$waterTemp)
md.Data$rainfall <- scale(md.Data$rainfall)
md.Data$roll_sum_24h <- scale(md.Data$roll_sum_24h)

md.Data$acce.netWeight <- scale(md.Data$acce.netWeight)

#### Remove data that is shorter than 1 day 
md.Data <- md.Data %>% 
  group_by(turtle.ID) %>% 
  arrange(time) %>% 
  split_at_gap(max_gap = 15, shortest_track = 1440) 

md.Data_iid <- split(md.Data, md.Data$turtle.ID)

md.Data_iid <- Filter(length, md.Data_iid)
md.Data_iid <- md.Data_iid[sapply(md.Data_iid, nrow) >= 288]

md.Data_iid %>% 
  bind_rows() %>% 
  group_by(ID) %>% 
  count() %>% 
  summarise(day = (n * 15)/(24*60)) %>% 
  summarise(mean = mean(day), 
            sd = sd(day), 
            max = max(day))

md.Data_iid <- lapply(md.Data_iid, function(x){
  
  if(nrow(x) > 0){
    
    x <- prepData(x, coord = NULL)
    
    x$ODBA.scaled <- log(x$ODBA + 0.001)

    return(x)
  }  
  
})

### Defining the model parameters ----
createCluster <- function(x){
  
  cluster <- mixmodCluster(x$ODBA.scaled, nbCluster = 2)
  
  hist(cluster, breaks = 100, main = NULL)
  title(main = x$turtle.ID[1])
  
  return(cluster)
  
}
emCluster <- lapply(md.Data_iid, createCluster)

createClusterTbl <- function(x){
  
  lapply(seq_along(x), function(i){
    
    emClusterMean <- c(x[[i]]@results[[1]]@parameters@mean)
    emClusterVar <- unlist(x[[i]]@results[[1]]@parameters@variance)
    emClusterSd <- sqrt(emClusterVar)
    emClusterTbl <- data.frame(ID = names(x)[i], 
                               mean = emClusterMean, 
                               var = emClusterVar,
                               sd = emClusterSd) %>% arrange(mean)
    emClusterTbl$state <- c(1:nrow(emClusterTbl))
    return(emClusterTbl)
    
  })
  
}
emClusterTbl <- createClusterTbl(emCluster) %>% 
  bind_rows() %>% 
  group_by(state) %>% 
  summarise(mean.min = min(mean), 
            mean.max = max(mean), 
            sd.min = min(sd), 
            sd.max = max(sd))

set.seed(1234)
seed <- sample(1:1000, length(md.Data_iid))

ODBAPar0 = list()
for(i in 1:length(seed)){
  set.seed(seed[i])
  n = 100
  ODBAPar0[[i]] <- data.frame(iter = c(1:n), 
                              mean1 = runif(n, min = emClusterTbl$mean.min[1], max = emClusterTbl$mean.max[1]), 
                              mean2 = runif(n, min = emClusterTbl$mean.min[2], max = emClusterTbl$mean.max[2]), 
                              sd1 = runif(n, min = emClusterTbl$sd.min[1], max = emClusterTbl$sd.max[1]), 
                              sd2 = runif(n, min = emClusterTbl$sd.min[2], max = emClusterTbl$sd.max[2]))
}

nbStates <- 2
stateNames <- c("stationary","mobile")

dist <- list(ODBA = "lnorm") 

#### Constrain ODBA parameters 
#### Mean(mobile) > Mean(stationary)
odbaDM0 <- matrix(c(1,0,0,0,0,0,
                    1,1,0,0,0,0,
                    0,0,1,0,0,0,
                    0,0,0,1,0,0,
                    0,0,0,0,1,0,
                    0,0,0,0,0,1),
                  nrow = 3*nbStates,byrow=TRUE,
                  dimnames=list(c(paste0("location_",1:nbStates),paste0("scale_",1:nbStates), paste0("zeromass_",1:nbStates)),
                                c("location_12:(Intercept)", "location_2", 
                                  paste0("scale_",1:nbStates,":(Intercept)"),
                                  paste0("zeromass_",1:nbStates,":(Intercept)"))))
odbaDM <- matrix(c(1,0,0,0,0,0,0,0,
                   1,1,"airTemp","roll_sum_24h",0,0,0,0,
                   0,0,0,0,1,0,0,0,
                   0,0,0,0,0,1,0,0,
                   0,0,0,0,0,0,1,0,
                   0,0,0,0,0,0,0,1),
                 nrow = 3*nbStates,byrow=TRUE,
                 dimnames=list(c(paste0("location_",1:nbStates),paste0("scale_",1:nbStates), paste0("zeromass_",1:nbStates)),
                               c("location_12:(Intercept)", "location_2", "location_2:airTemp", "location_2:roll_sum_24h",
                                 paste0("scale_",1:nbStates,":(Intercept)"),
                                 paste0("zeromass_",1:nbStates,":(Intercept)"))))
DM0 <- list(ODBA = odbaDM0)
DM <- list(ODBA = odbaDM)

#### Define the directions of the differences
odbaworkBounds0 <- matrix(c(-Inf, 0, rep(-Inf,4), 
                            rep(Inf,6)),
                          nrow = ncol(odbaDM0),
                          dimnames = list(colnames(odbaDM0), 
                                          c("lower", "upper")))

odbaworkBounds <- matrix(c(-Inf, 0, -10, -10, rep(-Inf,4),
                           rep(Inf,2), 10, 10, rep(Inf, 4)),
                         nrow = ncol(odbaDM),
                         dimnames = list(colnames(odbaDM),
                                         c("lower", "upper")))
workBounds0 <-list(ODBA = odbaworkBounds0)
workBounds <-list(ODBA = odbaworkBounds)

odbaBounds0 <- matrix(c(-10,0,
                        -10,0, 
                        0,Inf, 
                        0,Inf, 
                        0,1,
                        0,.01), 
                      nrow = ncol(odbaDM0), byrow = TRUE,
                      dimnames = list(colnames(odbaDM0), 
                                      c("lower", "upper")))
userBounds0 <- list(ODBA = odbaBounds0)

prior <- function(par){sum(dnorm(par,0,10000,log=TRUE))} ## prevents working parameters from straying along boundary

### Fitting models ----
hmmFits_ind_m1 = list()
hmmFits_ind_m2 = list()
for(i in 1:length(md.Data_iid)){
  
  hmmFits_ind_m1[[i]] = list()
  hmmFits_ind_m2[[i]] = list()
  
  names(hmmFits_ind_m1)[[i]] <- names(md.Data_iid)[[i]]
  names(hmmFits_ind_m2)[[i]] <- names(md.Data_iid)[[i]]
  
  zeromass1 <- nrow(md.Data_iid[[i]][md.Data_iid[[i]]$ODBA == 0,])/nrow(md.Data_iid[[i]]) 
  
  if(sum(md.Data_iid[[i]]$heavyRain=="yes")>10){
    
    formula <- ~ cosinor(hour, period = 24) + heavyRain
  
    }
  else{
    
    formula <- ~ cosinor(hour, period = 24)
    
  }
  
  for(j in 1:nrow(ODBAPar0[[i]])){
    
    tryCatch({

      Par <- c(as.numeric(ODBAPar0[[i]][j, 2:5]), zeromass1, 0) 
      Par <- list(ODBA = Par) ## dist == "lnorm" 

      hmmFits_ind_m1[[i]][[j]] <- fitHMM(md.Data_iid[[i]], 
                                         nbStates = nbStates, 
                                         dist = dist, 
                                         Par0 = Par, 
                                         formula = formula, 
                                         DM = DM0,
                                         workBounds = workBounds0, 
                                         userBounds = userBounds0, 
                                         prior = prior, 
                                         stateNames = stateNames)
      hmmFits_ind_m1[[i]][[j]]
      names(hmmFits_ind_m1[[i]])[[j]] <- paste0("mod", j)
      
      Par0_m2 <- getPar0(model = hmmFits_ind_m1[[i]][[j]],
                         formula = formula,
                         DM = DM, 
                         workBounds = workBounds, 
                         userBounds = userBounds0)

      hmmFits_ind_m2[[i]][[j]] <- fitHMM(md.Data_iid[[i]], nbStates = nbStates,
                                         dist = dist,
                                         Par0 = Par0_m2$Par, beta0 = Par0_m2$beta,
                                         formula = formula,
                                         DM = DM, 
                                         workBounds = workBounds,
                                         userBounds = userBounds0, 
                                         prior = prior, 
                                         stateNames = stateNames)
      hmmFits_ind_m2[[i]][[j]]
      names(hmmFits_ind_m2[[i]])[[j]] <- paste0("mod", j)
      
    }, error = function(e){cat("ERROR :",conditionMessage(e), "\n")})

  }
  
}

### Model selection ----
negLlk <- lapply(hmmFits_ind_m2, function(x){
  
  lapply(seq_along(x), function(i){
    
    x[[i]]$mod$minimum
    
  })

})

getBestMod <- function(x){
  
  mod = list()
  
  lapply(seq_along(x), function(i){

    modName <- names(x)[[i]]
    index <- as.integer(which.min(unlist(as.numeric(as.character(x[[i]])))))
    
    mod[[i]] <- hmmFits_ind_m2[[modName]][[index]]
    print(names(x))[[i]]

    return(mod[[i]])

  })    

}

bestMod <- getBestMod(negLlk)
ci <- lapply(bestMod, CIreal)

lapply(bestMod, plotPR) 

### Model evaluation ----
bestMod[[1]] 
bestMod[[1]]$CIbeta$ODBA
bestMod[[1]]$CIbeta$beta
plot(bestMod[[1]], plotCI = TRUE)

bestMod[[2]] 
bestMod[[2]]$CIbeta$ODBA
bestMod[[2]]$CIbeta$beta
plot(bestMod[[2]], plotCI = TRUE)

bestMod[[3]] 
bestMod[[3]]$CIbeta$ODBA
bestMod[[3]]$CIbeta$beta
plot(bestMod[[3]], plotCI = TRUE)

bestMod[[4]] 
bestMod[[4]]$CIbeta$ODBA
bestMod[[4]]$CIbeta$beta
plot(bestMod[[4]], plotCI = TRUE)

bestMod[[5]] 
bestMod[[5]]$CIbeta$ODBA
bestMod[[5]]$CIbeta$beta
plot(bestMod[[5]], plotCI = TRUE)

#### State processes ----
newData <- lapply(bestMod, function(x){
  
  newData <- x$data
  newData$state <- viterbi(x)
  
  return(newData)
  
})

lapply(newData, function(x){
  table(x$state)/nrow(x)
}) %>% 
  bind_rows() %>% 
  as.data.frame() %>% 
  dplyr::rename(c("s1" = "1", "s2" = "2")) %>% 
  summarise(mean = mean(s1), 
            sd = sd(s1))

#### Saving the hmm coded data 
# newData %>% 
#   bind_rows() %>% 
#   transform(time = as.character(format(time))) %>% 
#   write.csv("data/PLME_accelerometry_hmm.csv", row.names = FALSE)

#### State-dependent distribution and transition probabilities ----
odbaPara <- lapply(bestMod, function(x){
  
  id <- unique(x$data["turtle.ID"]) %>% t() %>% as.vector()
  
  state1_est <- x$CIreal$ODBA$est["location", "stationary"]  
  state1_lwr <- x$CIreal$ODBA$lower["location", "stationary"]
  state1_upr <- x$CIreal$ODBA$upper["location", "stationary"]
  
  state2_est <- x$CIreal$ODBA$est["location", "mobile"]
  state2_lwr <- x$CIreal$ODBA$lower["location", "mobile"]
  state2_upr <- x$CIreal$ODBA$upper["location", "mobile"]
  
  tp_11_est <- x$CIreal$gamma$est[1,1]
  tp_11_lwr <- x$CIreal$gamma$lower[1,1]
  tp_11_upr <- x$CIreal$gamma$upper[1,1]
  tp_12_est <- x$CIreal$gamma$est[1,2]
  tp_12_lwr <- x$CIreal$gamma$lower[1,2]
  tp_12_upr <- x$CIreal$gamma$upper[1,2]
  tp_21_est <- x$CIreal$gamma$est[2,1]
  tp_21_lwr <- x$CIreal$gamma$lower[2,1]
  tp_21_upr <- x$CIreal$gamma$upper[2,1]
  tp_22_est <- x$CIreal$gamma$est[2,2]
  tp_22_lwr <- x$CIreal$gamma$lower[2,2]
  tp_22_upr <- x$CIreal$gamma$upper[2,2]
  
  df <- data.frame(id, state1_est, state1_lwr, state1_upr, state2_est, state2_lwr, state2_upr, 
                   tp_11_est, tp_11_lwr, tp_11_upr, tp_12_est, tp_12_lwr, tp_12_upr, 
                   tp_21_est, tp_21_lwr, tp_21_upr, tp_22_est, tp_22_lwr, tp_22_upr)
  return(df)
  
})

#### Plotting -- State-dependent distribution ----
odbaParaData <- odbaPara %>% bind_rows()
odbaParaDataS1 <- cbind(odbaParaData[,1:4], state = "1")
names(odbaParaDataS1) <- c("id", "est", "lwr", "upr", "state")
odbaParaDataS2 <- cbind(odbaParaData[,c(1,5:7)], state = "2")
names(odbaParaDataS2) <- c("id", "est", "lwr", "upr", "state")
odbaParaData <- rbind(odbaParaDataS1, odbaParaDataS2)

#### Fig S1
odbaParaData %>% group_by(state) %>% 
  summarise(id = "Ensemble\nmean", 
            est = mean(est), 
            lwr = mean(lwr), 
            upr = mean(upr)) %>% 
  ggplot(aes(x = id, y = exp(est))) + 
  geom_point(cex = 5, col = "red") + 
  geom_errorbar(aes(ymin = exp(lwr), ymax = exp(upr)), width = 0, col = "red") + 
  geom_pointrange(data = odbaParaData, 
                  aes(x = id, y = exp(est), ymin = exp(lwr), ymax = exp(upr)), 
                  cex = 1, lty = 2) + 
  facet_grid(state~., scales = "free", 
             labeller = labeller(state = c("1" = "State 1", 
                                           "2" = "State 2"))) +
  theme_classic() + 
  theme(axis.title.x = element_blank(), 
        axis.ticks.x = element_blank()) + 
  ylab("ODBA (g)") + 
  theme(axis.title = element_text(size = 20)) +     
  theme(axis.text = element_text(size = 20)) +    
  theme(legend.text = element_text(size = 20)) + 
  theme(strip.text = element_text(size = 20))

#### Plotting -- State transition probabilities ----
tpBetas <- lapply(bestMod, function(x){
  tpBetas <- x$mle$beta
  tpBetas <- as.data.frame(t(tpBetas)) %>% 
    dplyr::select(Intercept.coef = '(Intercept)', 
                  Cos.coef = 'cosinorCos(hour, period = 24)', 
                  Sin.coef = 'cosinorSin(hour, period = 24)') %>% 
    mutate(Start_State = c(1,2)) %>% 
    mutate(End_State = c(2,1)) %>% 
    mutate(ID = unique(x[["data"]][["turtle.ID"]]))
})
tpBetas <- bind_rows(tpBetas)
tpLCL <- lapply(bestMod, function(x){
  tpLCL <- x$CIbeta$beta$lower
  tpLCL <- as.data.frame(t(tpLCL)) %>% 
    dplyr::select(Intercept.coef = '(Intercept)', 
                  Cos.coef = 'cosinorCos(hour, period = 24)', 
                  Sin.coef = 'cosinorSin(hour, period = 24)') %>% 
    mutate(Start_State = c(1,2)) %>% 
    mutate(End_State = c(2,1)) %>% 
    mutate(ID = unique(x[["data"]][["turtle.ID"]]))
})
tpLCL <- bind_rows(tpLCL)
tpUCL <- lapply(bestMod, function(x){
  tpUCL <- x$CIbeta$beta$upper
  tpUCL <- as.data.frame(t(tpUCL)) %>% 
    dplyr::select(Intercept.coef = '(Intercept)', 
                  Cos.coef = 'cosinorCos(hour, period = 24)', 
                  Sin.coef = 'cosinorSin(hour, period = 24)') %>% 
    mutate(Start_State = c(1,2)) %>% 
    mutate(End_State = c(2,1)) %>% 
    mutate(ID = unique(x[["data"]][["turtle.ID"]]))
})
tpUCL <- bind_rows(tpUCL)
tpData <- expand.grid(Hour.cov = seq(0,24,.1), 
                      Start_State = 1:2, 
                      End_State = 1:2)
tpData <- tpData %>% filter(Start_State != End_State)
tpFull <- merge(tpData, tpBetas, by = c("Start_State", "End_State"), all.x = TRUE) %>% 
  mutate(est = Intercept.coef + Cos.coef*cos(2*pi*Hour.cov/24) + Sin.coef*sin(2*pi*Hour.cov/24)) %>%  
  mutate(est = exp(est))
tpFull.LCL <- merge(tpData, tpLCL, by = c("Start_State", "End_State"), all.x = TRUE) %>% 
  mutate(est = Intercept.coef + Cos.coef*cos(2*pi*Hour.cov/24) + Sin.coef*sin(2*pi*Hour.cov/24)) %>%  
  mutate(est = exp(est))
tpFull.UCL <- merge(tpData, tpUCL, by = c("Start_State", "End_State"), all.x = TRUE) %>% 
  mutate(est = Intercept.coef + Cos.coef*cos(2*pi*Hour.cov/24) + Sin.coef*sin(2*pi*Hour.cov/24)) %>%  
  mutate(est = exp(est))
tp1 <- tpFull %>% 
  filter(Start_State == 1) %>% 
  dplyr::select(ID, Start_State, End_State, Hour.cov, est) %>% 
  pivot_wider(names_from = c(Start_State, End_State), values_from = est) %>% 
  rename(exp = '1_2') %>% 
  mutate(tp_2 = exp/(1+exp)) %>% 
  mutate(tp_1 = 1 - tp_2)
tp1$Start_State <- 1
tp2 <- tpFull %>% 
  filter(Start_State == 2) %>% 
  dplyr::select(ID, Start_State, End_State, Hour.cov, est) %>% 
  pivot_wider(names_from = c(Start_State, End_State), values_from = est) %>% 
  rename(exp = '2_1') %>% 
  mutate(tp_1 = exp/(1+exp)) %>% 
  mutate(tp_2 = 1 - tp_1)
tp2$Start_State <- 2
tp12 <- rbind(tp1, tp2)
tp1.lcl <- tpFull.LCL %>% 
  filter(Start_State == 1) %>% 
  dplyr::select(ID, Start_State, End_State, Hour.cov, est) %>% 
  pivot_wider(names_from = c(Start_State, End_State), values_from = est) %>% 
  rename(exp = '1_2') %>% 
  mutate(tp_2 = exp/(1+exp)) %>% 
  mutate(tp_1 = 1 - tp_2)
tp1.lcl$Start_State <- 1
tp2.lcl <- tpFull.LCL %>% 
  filter(Start_State == 2) %>% 
  dplyr::select(ID, Start_State, End_State, Hour.cov, est) %>% 
  pivot_wider(names_from = c(Start_State, End_State), values_from = est) %>% 
  rename(exp = '2_1') %>% 
  mutate(tp_1 = exp/(1+exp)) %>% 
  mutate(tp_2 = 1 - tp_1)
tp2.lcl$Start_State <- 2
tp12.lcl <- rbind(tp1.lcl, tp2.lcl)
tp1.ucl <- tpFull.UCL %>% 
  filter(Start_State == 1) %>% 
  dplyr::select(ID, Start_State, End_State, Hour.cov, est) %>% 
  pivot_wider(names_from = c(Start_State, End_State), values_from = est) %>% 
  rename(exp = '1_2') %>% 
  mutate(tp_2 = exp/(1+exp)) %>% 
  mutate(tp_1 = 1 - tp_2)
tp1.ucl$Start_State <- 1
tp2.ucl <- tpFull.UCL %>% 
  filter(Start_State == 2) %>% 
  dplyr::select(ID, Start_State, End_State, Hour.cov, est) %>% 
  pivot_wider(names_from = c(Start_State, End_State), values_from = est) %>% 
  rename(exp = '2_1') %>% 
  mutate(tp_1 = exp/(1+exp)) %>% 
  mutate(tp_2 = 1 - tp_1)
tp2.ucl$Start_State <- 2
tp12.ucl <- rbind(tp1.ucl, tp2.ucl)
tp12.cl <- merge(tp12.lcl, tp12.ucl, by = c("ID", "Start_State", "Hour.cov"))

#### Fig S2
ggplot(tp12, aes(x = Hour.cov)) +
  geom_line(aes(y = tp_1, color = "#1cade4"), size = 1) +
  geom_line(aes(y = tp_2, color = "#f1a806"), size = 1) +
  geom_ribbon(data = tp12.cl, aes(ymin = tp_1.x, ymax = tp_1.y),alpha = .3, fill = "#1cade4") +
  geom_ribbon(data = tp12.cl, aes(ymin = tp_2.x, ymax = tp_2.y),alpha = .3, fill = "#f1a806") +
  theme_classic() +
  scale_color_identity(name = "Behavioral\nState",
                       breaks = c("#1cade4", "#f1a806"),
                       labels = c("Stationary", "Mobile"),
                       guide = "legend")  +
  scale_linetype_identity(name = "Behavioral\nState",
                          breaks = c("longdash", "solid"),
                          labels = c("Stationary", "Mobile"),
                          guide = "legend") +
  xlab("Hour") +
  ylab("Transition probability") + 
  scale_x_continuous(breaks = (seq(0,24,4)),
                     labels = c(20, seq(0,20,4))) + 
  facet_grid(ID~Start_State, labeller = labeller(Start_State = c("1" = paste("P(Stationary ", sprintf("\u2192"), " X)"), 
                                                                 "2" = paste("P(Mobile ", sprintf("\u2192"), " X)")))) + 
  theme(axis.title = element_text(size = 12)) +     
  theme(axis.text = element_text(size = 12)) +    
  theme(legend.title = element_text(size = 12)) + 
  theme(legend.text = element_text(size = 12)) + 
  theme(strip.text = element_text(size = 12))

#### Sumamry statistics -- State-dependent distribution and transition probabilities ----
odbaParaSummary <- odbaPara %>% 
  bind_rows() %>% 
  summarise(state1_mean = exp(mean(state1_est)), 
            state1_sd = sd(exp(state1_est)), 
            state1_se = state1_sd/sqrt(length(state1_est)-1), 
            state1_lwr = exp(mean(state1_lwr)), 
            state1_upr = exp(mean(state1_upr)), 
            state2_mean = exp(mean(state2_est)),
            state2_sd = sd(exp(state2_est)), 
            state2_se = state2_sd/sqrt(length(state2_est)-1), 
            state2_lwr = exp(mean(state2_lwr)), 
            state2_upr = exp(mean(state2_upr)),
            tp11_mean = mean(tp_11_est), 
            tp11_sd = sd(tp_11_est), 
            tp11_se = tp11_sd/sqrt(length(tp_11_est)-1), 
            tp11_lwr = mean(tp_11_lwr), 
            tp11_upr = mean(tp_11_upr), 
            tp12_mean = mean(tp_12_est), 
            tp12_sd = sd(tp_12_est), 
            tp12_se = tp12_sd/sqrt(length(tp_12_est)-1), 
            tp12_lwr = mean(tp_12_lwr), 
            tp12_upr = mean(tp_12_upr), 
            tp21_mean = mean(tp_21_est), 
            tp21_sd = sd(tp_21_est), 
            tp21_se = tp21_sd/sqrt(length(tp_21_est)-1), 
            tp21_lwr = mean(tp_21_lwr), 
            tp21_upr = mean(tp_21_upr), 
            tp22_mean = mean(tp_22_est), 
            tp22_sd = sd(tp_22_est), 
            tp22_se = tp22_sd/sqrt(length(tp_22_est)-1), 
            tp22_lwr = mean(tp_22_lwr), 
            tp22_upr = mean(tp_22_upr))

#### Beta values ----
odbaBeta <- lapply(bestMod, function(x){
  
  airTemp <- x$CIbeta$ODBA$est[,"location_2:airTemp"]
  airTemp.lwr <- x$CIbeta$ODBA$lower[,"location_2:airTemp"]
  airTemp.upr <- x$CIbeta$ODBA$upper[,"location_2:airTemp"]
  
  rainfall <- x$CIbeta$ODBA$est[,"location_2:roll_sum_24h"]
  rainfall.lwr <- x$CIbeta$ODBA$lower[,"location_2:roll_sum_24h"]
  rainfall.upr <- x$CIbeta$ODBA$upper[,"location_2:roll_sum_24h"]
  
  if(c("heavyRainyes") %in% rownames(x$CIbeta$beta$est)){
    
    tp_heavyRain12 <- x$CIbeta$beta$est["heavyRainyes", "1 -> 2"]
    tp_heavyRain12.lwr <- x$CIbeta$beta$lower["heavyRainyes", "1 -> 2"]
    tp_heavyRain12.upr <- x$CIbeta$beta$upper["heavyRainyes", "1 -> 2"]
    
    tp_heavyRain21 <- x$CIbeta$beta$est["heavyRainyes", "2 -> 1"]
    tp_heavyRain21.lwr <- x$CIbeta$beta$lower["heavyRainyes", "2 -> 1"]
    tp_heavyRain21.upr <- x$CIbeta$beta$upper["heavyRainyes", "2 -> 1"]
    
  }
  else{ 
    
    tp_heavyRain12 <- tp_heavyRain21 <- NA
    tp_heavyRain12.lwr <- tp_heavyRain21.lwr <- NA
    tp_heavyRain12.upr <- tp_heavyRain21.upr <- NA
    
  }

  tp_hourSin12 <- x$CIbeta$beta$est["cosinorSin(hour, period = 24)", "1 -> 2"]
  tp_hourSin12.lwr <- x$CIbeta$beta$lower["cosinorSin(hour, period = 24)", "1 -> 2"]
  tp_hourSin12.upr <- x$CIbeta$beta$upper["cosinorSin(hour, period = 24)", "1 -> 2"]
  
  tp_hourSin21 <- x$CIbeta$beta$est["cosinorSin(hour, period = 24)", "2 -> 1"]
  tp_hourSin21.lwr <- x$CIbeta$beta$lower["cosinorSin(hour, period = 24)", "2 -> 1"]
  tp_hourSin21.upr <- x$CIbeta$beta$upper["cosinorSin(hour, period = 24)", "2 -> 1"]
  
  tp_hourCos12 <- x$CIbeta$beta$est["cosinorCos(hour, period = 24)", "1 -> 2"]
  tp_hourCos12.lwr <- x$CIbeta$beta$lower["cosinorCos(hour, period = 24)", "1 -> 2"]
  tp_hourCos12.upr <- x$CIbeta$beta$upper["cosinorCos(hour, period = 24)", "1 -> 2"]
  
  tp_hourCos21 <- x$CIbeta$beta$est["cosinorCos(hour, period = 24)", "2 -> 1"]
  tp_hourCos21.lwr <- x$CIbeta$beta$lower["cosinorCos(hour, period = 24)", "2 -> 1"]
  tp_hourCos21.upr <- x$CIbeta$beta$upper["cosinorCos(hour, period = 24)", "2 -> 1"]
  
  beta <- data.frame(airTemp, rainfall, 
                     tp_heavyRain12, tp_heavyRain21, 
                     tp_hourSin12, tp_hourSin21, tp_hourCos12, tp_hourCos21)
  betaCI <- data.frame(airTemp, airTemp.lwr, airTemp.upr, rainfall, rainfall.lwr, rainfall.upr, 
                       tp_heavyRain12, tp_heavyRain12.lwr, tp_heavyRain12.upr, 
                       tp_heavyRain21, tp_heavyRain21.lwr, tp_heavyRain21.upr,
                       tp_hourSin12, tp_hourSin12.lwr, tp_hourSin12.upr, 
                       tp_hourSin21, tp_hourSin21.lwr, tp_hourSin21.upr, 
                       tp_hourCos12, tp_hourCos12.lwr, tp_hourCos12.upr, 
                       tp_hourCos21, tp_hourCos21.lwr, tp_hourCos21.upr)
  
  return(list(beta, betaCI))
  
})
odbaBetaCI <- lapply(odbaBeta, function(x){x[[2]]}) %>% bind_rows()
odbaBeta <- lapply(odbaBeta, function(x){x[[1]]}) %>% bind_rows()

odbaBeta_mean <- sapply(odbaBeta, function(x) mean(x, na.rm = TRUE))
odbaBeta_sd <- sapply(odbaBeta, function(x) sd(x, na.rm = TRUE))

for(i in seq(1, ncol(odbaBetaCI), 3)){
  odbaBetaCI$newCol <- ifelse(odbaBetaCI[i+1]*odbaBetaCI[i+2]>0, 1, 0)
  names(odbaBetaCI)[names(odbaBetaCI) == "newCol"] <- paste0(names(odbaBetaCI)[i], ".sig")
}

tp_Beta <- odbaBetaCI %>%  
  mutate(index = rownames(.)) %>% 
  dplyr::select(c(index, starts_with("tp"))) %>% 
  pivot_longer(starts_with("tp"))
tp_Beta$trans <- gsub("\\D", "", tp_Beta$name)
tp_Beta$term <- gsub("(.+)_([[:alpha:]]+)([0-9].+)", "\\2", tp_Beta$name)
tp_Beta$term[tp_Beta$term == "heavyRain"] <- "heavyrain"
tp_Beta$est <- gsub("(.+)_([[:alpha:]]+)([0-9]+).(.+)", "\\4", tp_Beta$name)
tp_Beta$est[!tp_Beta$est %in% c("lwr", "upr", "sig")] <- "est"
tp_Beta$var <- gsub("(.+)([A-Z].+)","\\1", tp_Beta$term)
tp_Beta$fx <- gsub("(.+)([A-Z].+)","\\2", tp_Beta$term)

tp_Beta_hour <- tp_Beta %>% 
  filter(var == "hour") %>% 
  dplyr::select(-c(name, term)) %>% 
  pivot_wider(id_cols = c(index, trans, var, fx), 
              names_from = est, 
              values_from = value)
tp_Beta_rain <- tp_Beta %>% 
  filter(var == "heavyrain") %>% 
  dplyr::select(-c(name, term, fx)) %>% 
  pivot_wider(id_cols = c(index, trans, var), 
              names_from = est, 
              values_from = value)

#### Summary statistics -- Beta values ---- 
BetaSummary <- rbind(odbaBeta_mean, odbaBeta_sd)
BetaSummary <- t(BetaSummary)
BetaSummary <- as.data.frame(BetaSummary)
BetaSummary <- data.frame(term = rownames(BetaSummary), 
                          mean = BetaSummary$odbaBeta_mean, 
                          sd = BetaSummary$odbaBeta_sd)
tp_BetaSummary <- BetaSummary %>% filter(grepl("tp",term))
tp_BetaSummary$trans <- gsub("\\D", "", tp_BetaSummary$term)
tp_BetaSummary$term <- gsub("(.+)_([[:alpha:]]+)([0-9]+)", "\\2", tp_BetaSummary$term)
tp_BetaSummary$var <- gsub("(.+)([A-Z].+)","\\1", tp_BetaSummary$term)
tp_BetaSummary$fx <- gsub("(.+)([A-Z].+)","\\2", tp_BetaSummary$term)

s2_BetaSummary <- BetaSummary %>% filter(!grepl("tp",term))
s2_BetaSummary$fx <- gsub("(.+_.+)([A-Z].+)", "\\2", s2_BetaSummary$term)
s2_BetaSummary$term <- gsub("(.+_.+)([A-Z].+)", "\\1", s2_BetaSummary$term)

#### Plotting -- Beta values for emergent probabilities ----
#### Fig 2
odbaBetaPlot_1 <- s2_BetaSummary %>% 
  filter(term %in% c("airTemp")) %>% 
  ggplot(., aes(x = term, y = mean)) +  
  geom_point(aes(x = term, y = mean), cex = 5) + 
  geom_errorbar(aes(ymin = mean-sd, ymax = mean+sd), width = .1) + 
  geom_pointrange(data = cbind(term = "airTemp", odbaBetaCI), aes(y = airTemp, ymin = airTemp.lwr, ymax = airTemp.upr, col = factor(airTemp.sig, levels = c(0, 1))), 
                  position = sdamr::position_jitternudge(jitter.width = .2, nudge.x = .2, seed = 1101), 
                  cex = .8, linetype = 'dotted', alpha = .8, show.legend = TRUE) +
  scale_color_manual(values = c("1" = "darkblue", "0" = "lightblue"), 
                     labels = c("1" = "Significant effect", 
                                "0" = "Non-significant effect"), 
                     drop = FALSE) + 
  scale_y_continuous(position = "right",
                     limits = c(-(ceiling(max(abs(c(odbaBetaCI$airTemp.lwr, odbaBetaCI$airTemp.upr))))),
                                (ceiling(max(abs(c(odbaBetaCI$airTemp.lwr, odbaBetaCI$airTemp.upr))))))) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  facet_grid(~term, scales = "free_y", switch = "y", 
             labeller = labeller(term = c("airTemp" = "Air Temperature"))) + 
  theme_bw() + 
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank(), 
        legend.position = "bottom") + 
  theme(text = element_text(size = 12), 
        axis.text = element_text(size = 12), 
        strip.text = element_text(size = 12),
        legend.text = element_text(size = 12), 
        legend.title = element_blank())
odbaBetaPlot_2 <- s2_BetaSummary %>% 
  filter(term %in% c("rainfall")) %>% 
  ggplot(., aes(x = term, y = mean)) +  
  geom_point(aes(x = term, y = mean), cex = 5) + 
  geom_errorbar(aes(ymin = mean-sd, ymax = mean+sd), width = .1) + 
  geom_pointrange(data = cbind(term = "rainfall", odbaBetaCI), aes(y = rainfall, ymin = rainfall.lwr, ymax = rainfall.upr, col = factor(rainfall.sig)), 
                  position = sdamr::position_jitternudge(jitter.width = .2, nudge.x = .2, seed = 1101), 
                  cex = .8, linetype = 'dotted', alpha = .8) + 
  scale_color_manual(values = c("1" = "darkblue", "0" = "lightblue"), 
                     labels = c("1" = "Significant effect", 
                                "0" = "Non-significant effect")) + 
  scale_y_continuous(position = "right", 
                     limits = c(-(ceiling(max(abs(c(odbaBetaCI$rainfall.lwr, odbaBetaCI$rainfall.upr))))),
                                (ceiling(max(abs(c(odbaBetaCI$rainfall.lwr, odbaBetaCI$rainfall.upr))))))) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") + 
  facet_grid(~term, scales = "free_y", switch = "y", 
             labeller = labeller(term = c("rainfall" = "Rainfall"))) + 
  theme_bw() + 
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank(), 
        legend.position = "bottom") + 
  theme(text = element_text(size = 12), 
        axis.text = element_text(size = 12), 
        strip.text = element_text(size = 12),
        legend.text = element_text(size = 12), 
        legend.title = element_blank())
ggarrange(odbaBetaPlot_1, odbaBetaPlot_2, ncol=2, nrow=1, common.legend = TRUE, legend="bottom")  

#### Plotting -- Beta values for transition probabilities ---- 
#### Fig 3
tp_rain <- tp_BetaSummary %>% 
  filter(term == "heavyRain") %>% 
  ggplot(., aes(x = term, y = mean)) + 
  geom_pointrange(data = cbind(term = "heavyRain", tp_Beta_rain), 
                  aes(y = est, ymin = lwr, ymax = upr, col = factor(sig)), 
                  position = sdamr::position_jitternudge(jitter.width = .2, nudge.x = .2, seed = 1101), 
                  cex = .8, linetype = 'dotted', alpha = .8) + 
  scale_color_manual(values = c("1" = "darkblue", "0" = "lightblue"), 
                     labels = c("1" = "Significant effect", 
                                "0" = "Non-significant effect"), 
                     na.translate = FALSE, guide = "none") + 
  geom_point(aes(x = term, y = mean), cex = 5) + 
  geom_errorbar(aes(ymin = mean-sd, ymax = mean+sd), width = .1) + 
  scale_y_continuous(position = "right", 
                     limits = c(-(ceiling(max(abs(c(tp_Beta_rain$upr, tp_Beta_rain$lwr))))),
                                (ceiling(max(abs(c(tp_Beta_rain$upr, tp_Beta_rain$lwr))))))) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") + 
  facet_grid(term~trans, labeller = labeller(trans = c("12" = "stationary -> mobile", 
                                                       "21" = "mobile -> stationary"), 
                                             term = c("heavyRain" = "Heavy Rain")), 
             scales = "free_y", switch = "y") + 
  theme_bw() + 
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(), 
        axis.title.y = element_blank(), 
        legend.position = "none") + 
  theme(text = element_text(size = 12), 
        axis.text = element_text(size = 12), 
        strip.text = element_text(size = 12),
        legend.text = element_text(size = 12), 
        legend.title = element_blank())
tp_rain <- ggplotGrob(tp_rain)

tp_hour <- tp_BetaSummary %>% 
  filter(var == "hour") %>% 
  ggplot(., aes(x = var, y = mean, group = fx)) + 
  geom_pointrange(data = tp_Beta_hour, 
                  aes(y = est, ymin = lwr, ymax = upr, col = factor(sig), shape = fx), 
                  position = position_jitterdodge(jitter.width = .2, dodge.width = 1, seed = 1101), 
                  cex = .8, linetype = 'dotted', alpha = .8) + 
  scale_color_manual(values = c("1" = "darkblue", "0" = "lightblue"), 
                     labels = c("1" = "Significant effect", 
                                "0" = "Non-significant effect"), 
                     na.translate = FALSE) + 
  geom_point(aes(x = var, y = mean, shape = fx), 
             cex = 5, position = position_dodge(width = .5)) + 
  geom_errorbar(aes(ymin = mean-sd, ymax = mean+sd), 
                width = .1, position = position_dodge(width = .5)) + 
  scale_y_continuous(position = "right", 
                     limits = c(-(ceiling(max(abs(c(tp_Beta_hour$upr, tp_Beta_hour$lwr))))),
                                (ceiling(max(abs(c(tp_Beta_hour$upr, tp_Beta_hour$lwr))))))) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") + 
  facet_grid(var~trans, labeller = labeller(trans = c("12" = "stationary -> mobile", 
                                                      "21" = "mobile -> stationary"), 
                                            var = c("hour" = "Hour")),  
             scales = "free_y", switch = "y") + 
  scale_shape_manual(values = c(16,15)) +
  theme_bw() + 
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank(), 
        legend.title = element_blank(), 
        legend.position = "bottom") + 
  theme(text = element_text(size = 12), 
        axis.text = element_text(size = 12), 
        strip.text = element_text(size = 12),
        legend.text = element_text(size = 12), 
        legend.title = element_blank()) + 
  scale_shape_discrete(labels = c("Cos(x)", "Sin(x)"))
tp_hour <- ggplotGrob(tp_hour)

maxHeight = grid::unit.pmax(tp_hour$heights, tp_rain$heights)
tp_hour$heights <- tp_rain$heights <- as.list(maxHeight)
gridExtra::grid.arrange(tp_rain, tp_hour, nrow=2)

