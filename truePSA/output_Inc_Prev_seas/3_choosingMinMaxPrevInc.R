# sum(result[,1]<10&result[,2]<4)
# 
# sum(result[,2]>=2&result[,2]<4.5) #799
# sum(result[,1]>=1.2&result[,1]<10) #1857
# sum(result[,1]>=3&result[,1]<7) #978
# 
# sum((result[,2]>=2&result[,2]<4.5)&(result[,1]>=3&result[,1]<7))
# 
#setwd("/Users/sai/OneDrive/MORU/Projects/PSA of Savannakhet model/truePSA/output_Inc_Prev_seas") #mac

###finding out which values to use####
######################################
setwd("D:/OneDrive/MORU/Projects/PSA of Savannakhet model/truePSA/output_Inc_Prev_seas/") #windows
result <- read.csv('result/psa_run_2017-07-25 18_26_18.csv')[,-1]
# 
library(manipulate)
manipulate(
  #a <- sum((result[,2]>=prevmin&result[,2]<prevmax)&(result[,1]>=incmin&result[,1]<incmax))
  plot(sum((result[,2]>=prevmin&result[,2]<prevmax)&(result[,1]>=incmin&result[,1]<incmax)), main=paste(sum((result[,2]>=prevmin&result[,2]<prevmax)&(result[,1]>=incmin&result[,1]<incmax)))),
  prevmin=slider(2.1,4),
  prevmax=slider(4.1,5),
  incmin=slider(1,2.5),
  incmax=slider(2.5,3)
)

#chosen values
#prev 3.8 to 4.5
#inc 2.3 to 2.8
params_generated <- read.csv('result/psa_run_2017-07-25 18_26_18.csv')[,-1]
params_plausible_base <- params_generated[(params_generated[,2]>=3.48&params_generated[,2]<5)&(params_generated[,1]>=2.08&params_generated[,1]<2.814),]
write.csv(params_plausible_base,paste('result/plausibleParametersBase_',nrow(params_plausible_base),'.csv', sep = ''), row.names = FALSE)


#COMBINING FIXED Intervention parameters

#1. getting default intervention parameters
defaultParams <- c(
  API = 10,
  bh_max = 20,
  eta = 30,
  covEDAT0 = 25,
  covITN0 = 70,
  effITN = 30,
  covIRS0 = 0,
  effIRS = 15,
  muC = 1,
  muA = 1,
  muU = 1,
  percfail2018 = 5,
  percfail2019 = 15,
  percfail2020 = 30,
  EDATscale = 1,
  covEDATi = 70,
  ITNscale = 1,
  covITNi = 90,
  IRSscale = 1,
  covIRSi = 90,
  lossd = 30,
  dm = 6,
  cmda_1 = 50,
  cmda_2 = 50,
  cmda_3 = 50,
  tm_1 = 9,
  tm_2 = 10,
  tm_3 = 11,
  effv_1 = 75,
  effv_2 = 80,
  effv_3 = 92,
  vh = 90,
  MSATscale = 1,
  covMSATi = 90,
  MSATsensC = 99,
  MSATsensA = 87,
  MSATsensU = 4 
)

length(defaultParams) #37
dim(params_plausible_base) #151 39 #additional 2 columns are incidence & prevalence values

#2. combining default intervention parameters with generated plausible parameters for baseline
#that is Baseline+Fixed Intervention parameters
#baseline 1:14, intervention 15:37
baseline <- params_plausible_base[,-c(1,2,17:39)]
intervention <- data.frame(matrix(rep(defaultParams[15:37],nrow(baseline)),nrow(baseline), byrow = TRUE))
names(intervention) <- names(defaultParams)[15:37]

params_plausible_fixedInt <- rbind(cbind(baseline,intervention)[-1,],defaultParams) #adding the default parameters at 102nd row
dim(params_plausible_fixedInt) #102 37
write.csv(params_plausible_fixedInt, 'result/plausibleParamFinal.csv', row.names = FALSE)


